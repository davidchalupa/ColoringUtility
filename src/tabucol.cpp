#include <math.h>

#include "tabucol.h"

long long tabucol(graph G, refer colors, int alpha, int A, int B, long long t_max, int stage, refer *result, std::chrono::high_resolution_clock::time_point &start_time, long long time_limit)
{
    long long t = 0;
    long long cycles = 0;
    long long fitness, best_fitness, previous_best_fitness = 0;
    long long old_fitness;
    refer i = 0;
    refer best[MAX_VERTICES];

    bool aux_output = true;

    bool found = false;
    bool changed = false;

    tabu_base *tabu_searcher = new tabu_base();
    tabu_searcher->init(G, colors, alpha, A, B);
    tabu_searcher->generate_random();
    fitness = tabu_searcher->get_fitness();
    old_fitness = fitness;
    best_fitness = fitness;


//    printf("Colors: %u\n", colors);
//    printf("Current fitness: %lld", fitness);

    while (1)
    {
        if (t > t_max || std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count() > time_limit)
        {
//            if (aux_output) printf("\n");
            break;
        }

//        if (t > 0 && t % t_max == 0)
//        {
//            if (changed)
//            {
//                printf("\nBest fitness: %lld (changed)\n", best_fitness);
//            }
//            else
//            {
//                printf("\nBest fitness: %lld\n", best_fitness);
//            }
//            changed = false;
//            tabu_searcher->assign_solution(best);
//            if (best_fitness == previous_best_fitness)
//            {
//                cycles++;
//            }
//            else
//            {
//                cycles = 0;
//            }
//            previous_best_fitness = best_fitness;
//            t = 0;
//            t_max += ls_length;
//            stage = 2 - stage;
//        }

        fitness = tabu_searcher->get_fitness();

//        if (0 == t % 1000)
//        {
//            if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
//            {
//                putchar('\b');
//            }
//            if (aux_output) printf("%lld", fitness);
//            old_fitness = fitness;
//        }

        if (fitness >= G->m)
        {
//            if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
//            {
//                putchar('\b');
//            }
//            if (aux_output) printf("%lld\n", fitness);
//            printf("Found the following solution:\n");
//            for (i=0;i<G->n;i++)
//            {
//                fprintf(f, "%u", tabu_searcher->get_solution()[i]);
//                if (G->n-1 != i)
//                {
//                    fprintf(f, ",");
//                }
//            }
//            fprintf(f, "\n");
//            fprintf(f, "Iterations: %lld\n", t);
//            fprintf(f, "Final number of conflicts: %lld\n", G->m-tabu_searcher->compute_fitness(tabu_searcher->get_solution()));
            found = true;
            break;
        }

        if (stage == 0)
        {
            tabu_searcher->mutate_tabu();
        }
        else if (stage == 1)
        {
            tabu_searcher->mutate_tabu_rls_b();
        }
        else
        {
            tabu_searcher->mutate_tabu_rls();
        }


        if (best_fitness <= tabu_searcher->get_fitness())
        {
            best_fitness = tabu_searcher->get_fitness();
            for (i=0;i<G->n;i++)
            {
                best[i] = tabu_searcher->get_solution()[i];
                changed = true;
            }
        }

        t++;
    }

    // copying the result
    for (i=0;i<G->n;i++)
    {
        result[i] = tabu_searcher->get_solution()[i];
    }

    // this is to free the dynamically allocated memory
    tabu_searcher->deinit();

    delete(tabu_searcher);

    return (! found);
}
