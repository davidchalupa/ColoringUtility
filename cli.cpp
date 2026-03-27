#include "cli.h"

#include "compute.h"

#include <chrono>

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

int cli::start_cli(int argc, char **argv)
{
    sprintf(filename_output, "out.txt");

    srand((unsigned long long) time(NULL));

    long long time_limit = 600;

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
        else if (! strcmp(argv[param_index], "--out"))
        {
            if (param_index+1 >= argc)
            {
                printf("Invalid number of parameters.\n");
                return 1;
            }
            sprintf(filename_output, "%s", argv[param_index+1]);
            param_index++;
        }
        else if (! strcmp(argv[param_index], "--time-limit"))
        {
            if (param_index+1 >= argc)
            {
                printf("Invalid number of parameters.\n");
                return 1;
            }
            time_limit = atoll(argv[param_index+1]);
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

    auto start = std::chrono::high_resolution_clock::now();
    double total_time = 0;

    refer *coloring = new refer[G->n];

    compute(G, coloring, time_limit);

    auto end = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    printf("CPU time: %0.2lf s.\n", total_time);


    printf("Saving the coloring found to: %s...\n", filename_output);

    FILE *f;

    f = fopen(filename_output, "w");

    for (i=0;i<G->n;i++)
    {
        fprintf(f, "%u", coloring[i]);
        if (G->n-1 != i)
        {
            fprintf(f, ",");
        }
    }
    fprintf(f, "\n");

    fclose(f);

    delete[](coloring);

    if (NULL != G)
    {
        free_graph();
    }

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
    input_graph(source, "col");
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
