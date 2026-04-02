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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// helper to perform vector-matrix multiplication: y = A * x
void multiply_adj(graph G, double *x, double *y) {
    for (refer i = 0; i < G->n; i++)
    {
        y[i] = 0.0;
        for (refer j = 0; j < G->V[i].edgecount; j++)
        {
            refer neighbor = G->V[i].sibl[j];
            y[i] += x[neighbor];
        }
    }
}

double get_max_eigenvalue(graph G, int shift_mode, double shift_val, int iterations) {
    double *x = new double[G->n];
    double *y = new double[G->n];

    // initialize vector
    for (refer i = 0; i < G->n; i++)
    {
        x[i] = (double) rand() / RAND_MAX;
    }

    double lambda = 0;
    for (int it = 0; it < iterations; it++) {
        multiply_adj(G, x, y);

        // if shift_mode is active, we are calculating eigenvalues of (A - shift_val * I)
        if (shift_mode) {
            for (refer i = 0; i < G->n; i++)
            {
                y[i] -= shift_val * x[i];
            }
        }

        // calculate norm and normalize
        double norm = 0;
        for (refer i = 0; i < G->n; i++)
        {
            norm += y[i] * y[i];
        }
        norm = sqrt(norm);
        if (norm < 1e-9) break;

        for (refer i = 0; i < G->n; i++)
        {
            x[i] = y[i] / norm;
        }
        lambda = norm;
    }

    // Rayleigh quotient to get the actual sign for the shifted version
    double dot_top = 0, dot_bottom = 0;
    multiply_adj(G, x, y);
    if (shift_mode)
    {
        for (refer i = 0; i < G->n; i++) y[i] -= shift_val * x[i];
    }
    for (refer i = 0; i < G->n; i++)
    {
        dot_top += x[i] * y[i];
        dot_bottom += x[i] * x[i];
    }

    delete[](x);
    delete[](y);

    return dot_top / dot_bottom;
}

refer spectral_lower_bound(graph G) {
    if (G->n == 0)
    {
        return 0;
    }
    if (G->m == 0)
    {
        return 1;
    }

    // find lambda_max (usually around n * density)
    double l_max = get_max_eigenvalue(G, 0, 0, 30);


    // find lambda_min
    // we shift the matrix by l_max to make the most negative eigenvalue
    // the one with the largest absolute magnitude
    double shifted_max = get_max_eigenvalue(G, 1, l_max, 30);
    double l_min = shifted_max + l_max;

    //    printf("l_max: %0.3lf\n", l_max);
//    printf("l_min: %0.3lf\n", l_min);

    if (fabs(l_min) < 1e-6)
    {
        // safety check
        return (int) ceil(1.0 + l_max);
    }

    double bound = 1.0 + (l_max / fabs(l_min));

    return (refer) ceil(bound);
}

void compute(graph G, refer *coloring, long long time_limit)
{
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

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

    delete(algorithm_greedyclique_instance);

    // IG/RLS itself

    algorithm_IGRLS_instance = new algorithm_igcol(G, brelaz_colors);

    std::vector<unsigned long long> max_t_stag_configs = {100, 1000, 2500, 5000, 10000, 20000/*, 100000, 1000000*/};
    refer coloring_size;
    refer best_lower_bound = 0;
    for (const auto max_t_stag : max_t_stag_configs)
    {
        printf("Starting IG/RLS for graph coloring and maximum clique (max iter stag = %llu)...\n", max_t_stag);
        algorithm_IGRLS_instance->igcol(G, false, coloring, &clique_size, clique, greedy_clique_size, max_t_stag, start, time_limit);
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

    printf("Attempting to use Hoffman spectral lower bound...\n");

    refer new_lower_bound = spectral_lower_bound(G);

    if (new_lower_bound > best_lower_bound)
    {
        best_lower_bound = new_lower_bound;
        printf("Found a new lower bound of %d colors.\n", best_lower_bound);
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
        if (0 == tabucol(G, k, 0, 0, 1, t_max, 0, potential_coloring, start, time_limit))
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
        if (0 == tabucol(G, k, 6, 10, 1, t_max, 0, potential_coloring, start, time_limit))
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
        if (0 == tabucol(G, k, 6, 10, 1, t_max, 1, potential_coloring, start, time_limit))
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
//        if (0 == tabucol(G, k, 6, 10, 1, t_max, 2, coloring, start, time_limit))
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
