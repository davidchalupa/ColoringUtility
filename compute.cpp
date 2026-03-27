#include "compute.h"

#include "algorithm_brelaz.h"
#include "algorithm_greedyclique.h"
#include "algorithm_igcol.h"
#include "tabucol.h"
#include "h2col.h"

#include <vector>

#include <chrono>

refer count_colors(graph G, refer *result)
{
    refer v,max_label;
    max_label = 0;
    for (v=0;v<G->n;v++)
    {
        if (max_label < result[v])
        {
            max_label = result[v];
        }
    }
    return max_label;
}

void compute(graph G, refer *coloring, long long time_limit)
{
    auto start = std::chrono::high_resolution_clock::now();

    algorithm_brelaz *algorithm_brelaz_instance = new algorithm_brelaz();

    refer brelaz_colors = algorithm_brelaz_instance->brelaz_with_heap(G, coloring);

    printf("Found a solution with %d colors with Brelaz's heuristic.\n", brelaz_colors);

    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
    {
        printf("Time limit exceeded.\n");
        return;
    }

    refer best_coloring_size = brelaz_colors;

    int brelaz_extra_trials_no_improve = 5;

    printf("Attempting extra trials...\n");

    int trial_no_improve_index = 0;
    while (trial_no_improve_index < brelaz_extra_trials_no_improve)
    {
        brelaz_colors = algorithm_brelaz_instance->brelaz_with_heap(G, coloring);
        if (brelaz_colors < best_coloring_size)
        {
            best_coloring_size = brelaz_colors;
            printf("Found a solution with %d colors with Brelaz's heuristic.\n", brelaz_colors);
            trial_no_improve_index = 0;
        }
        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
        {
            printf("Time limit exceeded.\n");
            return;
        }
        trial_no_improve_index++;
    }

    delete(algorithm_brelaz_instance);

    algorithm *algorithm_IGRLS_instance;
    algorithm *algorithm_greedyclique_instance;

    refer greedy_clique_size;
    refer *clique;

    refer clique_size;

    refer v;

    srand((unsigned long long) time(NULL));

    clique = new refer[G->n];

    algorithm_greedyclique_instance = new algorithm_greedyclique();

    // initializing the clique with greedy algorithm
    greedy_clique_size = algorithm_greedyclique_instance->greedy_clique(G, clique);

    // IG/RLS itself

    algorithm_IGRLS_instance = new algorithm_igcol(G, brelaz_colors);

    std::vector<unsigned long long> max_t_stag_configs = {100, 1000, 2500, 5000, 10000, 20000/*, 100000, 1000000*/};
    refer coloring_size;
    refer best_lower_bound = 0;
    for (const auto max_t_stag : max_t_stag_configs)
    {
        printf("Starting IG/RLS for graph coloring and maximum clique (max iter stag = %llu)...\n", max_t_stag);
        algorithm_IGRLS_instance->igcol(G, false, coloring, &clique_size, clique, greedy_clique_size, max_t_stag);
        coloring_size = count_colors(G, coloring);

        bool should_continue = false;
        if (clique_size > best_lower_bound)
        {
            best_lower_bound = clique_size;
            printf("Found a new lower bound of %d colors.\n", best_lower_bound);
            should_continue = true;
        }
        if (coloring_size < best_coloring_size)
        {
            best_coloring_size = coloring_size;
            printf("Found a solution with %d colors with IG/RLS.\n", best_coloring_size);
            should_continue = true;
        }

        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
        {
            printf("Time limit exceeded.\n");
            return;
        }

        // rather forcing continuing with small numbers of iterations
        if (max_t_stag <= 2500)
        {
            should_continue = true;
        }

        if (best_lower_bound == best_coloring_size || ! should_continue)
        {
            break;
        }
    }

    delete(algorithm_IGRLS_instance);

    if (best_lower_bound == best_coloring_size)
    {
        printf("Found an optimal coloring with %d colors!\n", best_coloring_size);
        return;
    }

    long long t;
    int i;

    int tabu_tenure;
    unsigned long long t_max;
    unsigned long long ls_length;
    tabu_tenure = 0;
    t_max = 100000;
    ls_length = 10000;

    refer *potential_coloring = new refer[G->n];

    refer k = best_coloring_size - 1;

    // ToDo: or time limit
    while (best_lower_bound <= k)
    {
        printf("Starting k-fixed local search (k = %d)...\n", k);
        if (0 == tabucol(G, k, 0, 0, 1, t_max, 0, potential_coloring))
        {
            printf("Found a coloring with %d colors.\n", k);
            best_coloring_size = k;
            for (refer v = 0; v < G->n; v++)
            {
                coloring[v] = potential_coloring[v];
            }
            k--;
        }
        else
        {
            break;
        }

        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
        {
            printf("Time limit exceeded.\n");
            return;
        }
    }

    if (best_lower_bound == best_coloring_size)
    {
        printf("Found an optimal coloring with %d colors!\n", best_coloring_size);
        delete[](potential_coloring);
        return;
    }

    while (best_lower_bound <= k)
    {
        printf("Starting TabuCol with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
        if (0 == tabucol(G, k, 6, 10, 1, t_max, 0, potential_coloring))
        {
            printf("Found a coloring with %d colors.\n", k);
            for (refer v = 0; v < G->n; v++)
            {
                coloring[v] = potential_coloring[v];
            }
            best_coloring_size = k;
            k--;
            continue;
        }

        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
        {
            printf("Time limit exceeded.\n");
            return;
        }

        printf("Starting RLS_b with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
        if (0 == tabucol(G, k, 6, 10, 1, t_max, 1, potential_coloring))
        {
            printf("Found a coloring with %d colors.\n", k);
            for (refer v = 0; v < G->n; v++)
            {
                coloring[v] = potential_coloring[v];
            }
            best_coloring_size = k;
            k--;
            continue;
        }

        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() > time_limit)
        {
            printf("Time limit exceeded.\n");
            return;
        }

//        printf("Starting RLS_a with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
//        if (0 == tabucol(G, k, 6, 10, 1, t_max, 2, coloring))
//        {
//            printf("Found a coloring with %d colors.\n", k);
//            best_coloring_size = k;
//            k--;
//            continue;
//        }
        t_max *= 10;
    }

    if (best_lower_bound == best_coloring_size)
    {
        printf("Found an optimal coloring with %d colors!\n", best_coloring_size);
        delete[](potential_coloring);
        return;
    }
    else
    {
        printf("Found a coloring with %d colors.\n", best_coloring_size);
    }

    //if (0 == hyblscol(G, k, t_max, t_max, &t, f))
    // TabuCol with static tabu tenure
    //if (0 == tabucol(G, k, 0, 0, 1, t_max, t_max, &t, f))
    // H2Col with static tabu tenure
    //if (0 == h2col(G, k, 6, 10, 0, t_max, &t, f))
    // the new coloring algorithm
    //if (0 == newcol(G, k, 6, 10, tabu_tenure, ls_length, t_max, ils_cycles_max, improvement_cycles_max, &t, f))

    delete[](potential_coloring);
}
