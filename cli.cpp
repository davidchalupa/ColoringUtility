#include "cli.h"

#include "algorithm_brelaz.h"
#include "algorithm_greedyclique.h"
#include "algorithm_igcol.h"
#include "tabucol.h"
#include "h2col.h"

#include <vector>

#include <QTime>

cli::cli()
{
    G = NULL;
}

void cli::sleep(int milisec)
{
    // ToDo: non-Qt implementation
}

void cli::usage()
{
    printf("This is a solver for the graph coloring instances.\n");
    printf("Usage: ColoringUtility --in <input_file.col> --k <the_number_of_colors> --out <output_file.txt>\n\n");
    printf("       --in <input_file.col>: switch specifies the input file in COL format\n");
    printf("       --out <output_file.txt>: switch specifies the output file name\n");
}

void cli::process_params()
{
}

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

void cli::compute()
{
    algorithm_brelaz *algorithm_brelaz_instance = new algorithm_brelaz();
    refer *coloring = new refer[G->n];

    refer brelaz_colors = algorithm_brelaz_instance->brelaz_with_heap(G, coloring);

    printf("Found a solution with %d colors with Brelaz's heuristic.\n", brelaz_colors);

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

    std::vector<unsigned long long> max_t_stag_configs = {100, 1000, 2500, 10000/*, 100000, 1000000*/};
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

        if (best_lower_bound == best_coloring_size || ! should_continue)
        {
            break;
        }
    }

    delete(algorithm_IGRLS_instance);

    delete[](coloring);

    if (best_lower_bound == best_coloring_size)
    {
        printf("Found an optimal coloring with %d colors!\n", coloring_size);
        return;
    }

    long long t;
    int i;

    tabu_tenure = 0;
    t_max = 100000;
    ls_length = 10000;


    k = best_coloring_size - 1;

    // ToDo: or time limit
    while (best_lower_bound <= k)
    {
        printf("Starting k-fixed local search (k = %d)...\n", k);
        if (0 == tabucol(G, k, 0, 0, 1, t_max, t_max, 0, &t))
        {
            printf("Found a coloring with %d colors.\n", k);
            best_coloring_size = k;
            k--;
        }
        else
        {
            break;
        }
    }

    if (best_lower_bound == best_coloring_size)
    {
        printf("Found an optimal coloring with %d colors!\n", best_coloring_size);
        return;
    }

    while (best_lower_bound <= k)
    {
        printf("Starting TabuCol with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
        if (0 == tabucol(G, k, 6, 10, 1, ls_length, t_max, 0, &t))
        {
            printf("Found a coloring with %d colors.\n", k);
            best_coloring_size = k;
            k--;
            continue;
        }
        printf("Starting RLS_b with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
        if (0 == tabucol(G, k, 6, 10, 1, ls_length, t_max, 1, &t))
        {
            printf("Found a coloring with %d colors.\n", k);
            best_coloring_size = k;
            k--;
            continue;
        }
//        printf("Starting RLS_a with dynamic tabu tenure (k = %d, t_max = %lld)...\n", k, t_max);
//        if (0 == tabucol(G, k, 6, 10, 1, ls_length, t_max, 2, &t))
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
}


int cli::start_cli(int argc, char **argv)
{
    sprintf(filename_output, "out.txt");

    srand((unsigned long long) time(NULL));

    if (argc < 2)
    {
        usage();
        return 1;
    }

    int param_index = 1;
    while (param_index < argc)
    {
        if (! strcmp(argv[param_index], "--in"))
        {
            if (param_index+1 >= argc)
            {
                printf("Invalid number of parameters.\n");
                return 1;
            }
            sprintf(filename, "%s", argv[param_index+1]);
            param_index++;
        }
        else
        {
            printf("Unknown argument %s.\n", argv[param_index]);
            return 1;
        }

        param_index++;
    }

    choose_instance();

    compute_statistics();

    srand((unsigned) time(0));

    FILE *f;

    f = fopen(filename_output, "w");

    QTime qTime;
    qTime.start();
    double total_time = 0;

    compute();

    total_time += qTime.elapsed() / 1000.0;
    printf("CPU time: %0.2lf s.\n", total_time);

    if (NULL != G)
    {
        free_graph();
    }

    fclose(f);

    return 0;
}

int cli::generate_instance()
{
    if (get_graph() != NULL)
    {
        free_graph();
    }
    generate_graph_BA_model(w,n_max);
    G = get_graph();

    return 0;
}

int cli::generate_instance_udg()
{
    if (get_graph() != NULL)
    {
        free_graph();
    }
    generate_graph_UDG(n_max,range,grid);
    G = get_graph();

    return 0;
}

int cli::choose_instance()
{
    if ((source = fopen(filename,"r")) == NULL)
    {
        printf("Error: An error occured during opening of the file.");
        return 1;
    }
    if (get_graph() != NULL)
    {
        free_graph();
    }
    input_graph(source);
    fclose(source);
    G = get_graph();

    return 0;
}

int cli::compute_statistics()
{
//    refer i, maxdeg;

//    min_deg = statistics::min_degree(G);
//    max_deg = statistics::max_degree(G);
//    avg_deg = statistics::average_degree(G);
//    stdev_deg = statistics::degree_stdev(G);
//    triangles = statistics::triangles(G);

    printf("G: %s\n",filename);
    printf("|V| = %u, |E| = %lu, d = %0.3lf\n",G->n,G->m,(double)(2*G->m)/(double)(G->n)/(double)(G->n-1));
//    printf("min degree = %ld, max degree = %ld\n",min_deg,max_deg);
//    printf("average degree = %0.3lf, degree stdev = %0.3lf\n",avg_deg,stdev_deg);
//    printf("number of triangles = %lu\n",triangles);

//    printf("degree distribution: [");
//    maxdeg = statistics::degree_distribution(G, degree_distrib);
//    for (i=0;i<=maxdeg;i++)
//    {
//        printf("%u", degree_distrib[i]);
//        if (i != maxdeg)
//        {
//            putchar(',');
//        }
//    }
//    printf("]\n");

    return 0;
}
