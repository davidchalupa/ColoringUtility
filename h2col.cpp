/*
 * This one a pure H2Col to sample promising local optima.
 */

#include <math.h>

#include "time.h"
#include "common.h"
#include "h2col.h"
#include "random_generator.h"

#define MAX_POPULATION 3
#define IMPOSSIBLE MAX_VERTICES + 1

typedef tabu_base *p_tabu_base;

refer dist_mtx[MAX_COLORS][MAX_COLORS];
refer M[MAX_COLORS];
refer sigma[MAX_COLORS];
refer cP1[MAX_COLORS];
refer cP2[MAX_COLORS];

long partition_distance_heuristic(graph G, refer colors, refer *P1, refer *P2)
{
    long Sigma = 0;
    refer i,j,K,x;
    K = colors;

    // M, sigma, cP1 and cP2 to 0 vectors
    for (i=1;i<=K;i++)
    {
        M[i] = 0;
        sigma[i] = 0;
        cP1[i] = 0;
        cP2[i] = 0;
    }

    for (x=0;x<G->n;x++)
    {
        i = P1[x];
        j = P2[x];
        dist_mtx[i][j] = 0;
    }

    for (x=0;x<G->n;x++)
    {
        i = P1[x];
        j = P2[x];
        dist_mtx[i][j]++;
        cP1[i]++;
        cP2[j]++;
        if (dist_mtx[i][j] > M[i])
        {
            M[i] = dist_mtx[i][j];
            sigma[i] = j;
        }
    }

    for (i=1;i<=K;i++)
    {
        if (0 == M[i])
        {
            continue;
        }

        if (3*M[i] <= cP1[i] + cP2[sigma[i]])
        {
            return IMPOSSIBLE;
        }
        Sigma += dist_mtx[i][sigma[i]];
    }

    return Sigma;
}

bool are_equal(graph G, refer *P1, refer *P2)
{
    refer i;

    for (i=0;i<G->n;i++)
    {
        if (P1[i] != P2[i])
        {
            return false;
        }
    }

    return true;
}

long long h2col(graph G, refer colors, int alpha, int A, int B, long long ls_length, long long *output_t, FILE *f)
{
    long long t = 0, g = 0;
    long long j;
    long long max_fitness;
    long long min_fitness;
    long long elite1_fitness;
    long long elite2_fitness;
    long long best_fitness;
    long long iter_cycle = 10;
    long long final_conflicts;

    int max_fitness_index = 0;
    long population_size = 2;
    int i;
    refer similarity = 0;

    bool aux_output = true;

    bool found = false;

    p_tabu_base tabu_searchers[MAX_POPULATION];
    long long fitnesses[MAX_POPULATION];
    for (i=0;i<population_size;i++)
    {
        tabu_searchers[i] = new tabu_base();
        tabu_searchers[i]->init(G, colors, alpha, A, B);
        {
            tabu_searchers[i]->generate_random();
            fitnesses[i] = tabu_searchers[i]->get_fitness();
            for (j=0;j<ls_length;j++)
            {
                tabu_searchers[i]->mutate_tabu();
            }
            tabu_searchers[i]->assign_solution(tabu_searchers[i]->get_best_solution());
        }
        fitnesses[i] = tabu_searchers[i]->get_fitness();
    }

    // offspring allocation
    p_tabu_base offspring1;
    offspring1 = new tabu_base();
    offspring1->init(G, colors, alpha, A, B);

    p_tabu_base offspring2;
    offspring2 = new tabu_base();
    offspring2->init(G, colors, alpha, A, B);

    // elite generation
    p_tabu_base elite1;
    elite1 = new tabu_base();
    elite1->init(G, colors, alpha, A, B);
    elite1->generate_random();
    elite1_fitness = elite1->get_fitness();

    p_tabu_base elite2;
    elite2 = new tabu_base();
    elite2->init(G, colors, alpha, A, B);
    elite2->generate_random();
    elite2_fitness = elite2->get_fitness();

    p_tabu_base best;
    best = new tabu_base();
    best->init(G, colors, alpha, A, B);
    best->generate_random();
    best_fitness = best->get_fitness();

    while (1)
    {
        // fitnesses and their maximum
        max_fitness = 0;
        min_fitness = G->m;
        for (i=0;i<population_size;i++)
        {
            if (max_fitness < fitnesses[i])
            {
                max_fitness_index = i;
                max_fitness = fitnesses[i];
            }
            if (min_fitness > fitnesses[i])
            {
                min_fitness = fitnesses[i];
            }
        }

        // auxiliary output
        if (aux_output) printf("Generation %lld: %lld, %lld; %lld. LS length: %ld\n", g, fitnesses[0], fitnesses[1], best_fitness, ls_length);

        // optimum handling
        if (max_fitness >= G->m)
        {
            final_conflicts = G->m-tabu_searchers[max_fitness_index]->compute_fitness(tabu_searchers[max_fitness_index]->get_solution());
            fprintf(f, "Found the following solution:\n");
            for (i=0;i<G->n;i++)
            {
                fprintf(f, "%u", best->get_solution()[i]);
                if (G->n-1 != i)
                {
                    fprintf(f, ",");
                }
            }
            fprintf(f, "\n");
            fprintf(f, "Generations: %lld\n", g);
            fprintf(f, "Iterations: %lld\n", t);
            fprintf(f, "Final number of conflicts: %lld\n", G->m-best->compute_fitness(best->get_solution()));
            found = true;

            found = true;
            break;
        }

        // selection of parents and crossover
        offspring1->create_by_greedy_crossover(tabu_searchers[0]->get_solution(), tabu_searchers[1]->get_solution());
        offspring2->create_by_greedy_crossover(tabu_searchers[1]->get_solution(), tabu_searchers[0]->get_solution());

        // tabu search on the offspring
        for (j=0;j<ls_length;j++)
        {
            offspring1->mutate_tabu();
            t++;
        }
        // tabu search on the offspring
        for (j=0;j<ls_length;j++)
        {
            offspring2->mutate_tabu();
            t++;
        }

        if (elite1_fitness <= offspring1->get_best_solution_fitness())
        {
            elite1->assign_solution(offspring1->get_best_solution());
            elite1_fitness = elite1->get_fitness();
        }
        if (elite1_fitness <= offspring2->get_best_solution_fitness())
        {
            elite1->assign_solution(offspring2->get_best_solution());
            elite1_fitness = elite1->get_fitness();
        }

        tabu_searchers[0]->assign_solution(offspring1->get_best_solution());
        tabu_searchers[1]->assign_solution(offspring2->get_best_solution());

        fitnesses[0] = tabu_searchers[0]->get_fitness();
        fitnesses[1] = tabu_searchers[1]->get_fitness();

        if (best_fitness <= elite1_fitness)
        {
            best->assign_solution(elite1->get_solution());
            best_fitness = elite1_fitness;
        }

        if (t % iter_cycle == 0)
        {
            tabu_searchers[0]->assign_solution(elite2->get_solution());
            fitnesses[0] = tabu_searchers[0]->get_fitness();
            elite2->assign_solution(elite1->get_solution());
            elite1->generate_random();
        }
        for (i=0;i<population_size;i++)
        {
            fitnesses[i] = tabu_searchers[i]->get_fitness();
        }
        elite1_fitness = elite1->get_fitness();
        elite2_fitness = elite2->get_fitness();

        g++;
    }

    for (i=0;i<population_size;i++)
    {
        delete(tabu_searchers[i]);
    }
    delete(offspring1);
    delete(offspring2);
    delete(elite1);
    delete(elite2);
    delete(best);

    *output_t = t;

    return (! found);
}

long long newcol(graph G, refer colors, int alpha, int A, int B, long long ls_length, long short_ils_length, long ils_cycles_max, long improvement_cycles_max, long long *output_t, FILE *f)
{
    long long t = 0, g = 0;
    long long j;
    long long max_fitness;
    long long min_fitness;
    long long elite1_fitness;
    long long elite2_fitness;
    long long best_fitness;
    long long iter_cycle = 10;
    long long final_conflicts;

    long long short_ls_length = ls_length;

    int max_fitness_index = 0;
    long population_size = 2;
    int i;
    int stag;
    int stag_max = 1000;
    refer similarity = 0;

    bool aux_output = true;

    bool found = false;

    p_tabu_base tabu_searchers[MAX_POPULATION];
    long long fitnesses[MAX_POPULATION];
    for (i=0;i<population_size;i++)
    {
        tabu_searchers[i] = new tabu_base();
        tabu_searchers[i]->init(G, colors, alpha, A, 0);
        {
            tabu_searchers[i]->generate_random();
            fitnesses[i] = tabu_searchers[i]->get_fitness();
            for (j=0;j<short_ls_length;j++)
            {
                tabu_searchers[i]->mutate_tabu();
            }
            tabu_searchers[i]->assign_solution(tabu_searchers[i]->get_best_solution());
        }
        fitnesses[i] = tabu_searchers[i]->get_fitness();
    }

    // offspring allocation
    p_tabu_base offspring1;
    offspring1 = new tabu_base();
    offspring1->init(G, colors, alpha, A, 0);

    p_tabu_base offspring2;
    offspring2 = new tabu_base();
    offspring2->init(G, colors, alpha, A, 0);

    // elite generation
    p_tabu_base elite1;
    elite1 = new tabu_base();
    elite1->init(G, colors, alpha, A, 0);
    elite1->generate_random();
    elite1_fitness = elite1->get_fitness();

    p_tabu_base elite2;
    elite2 = new tabu_base();
    elite2->init(G, colors, alpha, A, 0);
    elite2->generate_random();
    elite2_fitness = elite2->get_fitness();

    p_tabu_base best;
    best = new tabu_base();
    best->init(G, colors, alpha, A, 0);
    best->generate_random();
    best_fitness = best->get_fitness();

    tabu_base *tabu_searcher = new tabu_base();
    tabu_searcher->init(G, colors, 0, 0, B);
    tabu_searcher->assign_solution(best->get_solution());
    long long fitness = tabu_searcher->get_fitness();
    long long old_fitness = fitness;
    long long previous_best_fitness;

    // PHASE 1 - short-term H2Col with dynamic tabu tenure

    stag = 0;
    while (1)
    {
        // fitnesses and their maximum
        max_fitness = 0;
        min_fitness = G->m;
        for (i=0;i<population_size;i++)
        {
            if (max_fitness < fitnesses[i])
            {
                max_fitness_index = i;
                max_fitness = fitnesses[i];
            }
            if (min_fitness > fitnesses[i])
            {
                min_fitness = fitnesses[i];
            }
        }

        // auxiliary output
        if (aux_output) printf("Generation %lld: %lld, %lld; %lld. LS length: %ld\n", g, fitnesses[0], fitnesses[1], best_fitness, short_ls_length);

        // optimum handling
        if (max_fitness >= G->m)
        {
            final_conflicts = G->m-tabu_searchers[max_fitness_index]->compute_fitness(tabu_searchers[max_fitness_index]->get_solution());
            fprintf(f, "Found the following solution:\n");
            for (i=0;i<G->n;i++)
            {
                fprintf(f, "%u", best->get_solution()[i]);
                if (G->n-1 != i)
                {
                    fprintf(f, ",");
                }
            }
            fprintf(f, "\n");
            fprintf(f, "Generations: %lld\n", g);
            fprintf(f, "Iterations: %lld\n", t);
            fprintf(f, "Final number of conflicts: %lld\n", G->m-best->compute_fitness(best->get_solution()));
            found = true;

            break;
        }

        similarity = partition_distance_heuristic(G, colors, tabu_searchers[0]->get_solution(), tabu_searchers[1]->get_solution());
        if (IMPOSSIBLE != similarity && similarity >= 99*G->n/100 || stag >= stag_max)
        {
            final_conflicts = G->m-tabu_searchers[max_fitness_index]->compute_fitness(tabu_searchers[max_fitness_index]->get_solution());
            break;
        }

        // selection of parents and crossover
        offspring1->create_by_greedy_crossover(tabu_searchers[0]->get_solution(), tabu_searchers[1]->get_solution());
        offspring2->create_by_greedy_crossover(tabu_searchers[1]->get_solution(), tabu_searchers[0]->get_solution());

        // tabu search on the offspring
        for (j=0;j<short_ls_length;j++)
        {
            offspring1->mutate_tabu();
            t++;
        }
        // tabu search on the offspring
        for (j=0;j<short_ls_length;j++)
        {
            offspring2->mutate_tabu();
            t++;
        }

        if (elite1_fitness <= offspring1->get_best_solution_fitness())
        {
            elite1->assign_solution(offspring1->get_best_solution());
            elite1_fitness = elite1->get_fitness();
        }
        if (elite1_fitness <= offspring2->get_best_solution_fitness())
        {
            elite1->assign_solution(offspring2->get_best_solution());
            elite1_fitness = elite1->get_fitness();
        }

        tabu_searchers[0]->assign_solution(offspring1->get_best_solution());
        tabu_searchers[1]->assign_solution(offspring2->get_best_solution());

        fitnesses[0] = tabu_searchers[0]->get_fitness();
        fitnesses[1] = tabu_searchers[1]->get_fitness();

        if (best_fitness <= elite1_fitness)
        {
            best->assign_solution(elite1->get_solution());
            best_fitness = elite1_fitness;
            stag = 0;
        }
        else
        {
            stag++;
        }

        if (g > 0 && g % iter_cycle == 0)
        {
            tabu_searchers[0]->assign_solution(elite2->get_solution());
            fitnesses[0] = tabu_searchers[0]->get_fitness();
            elite2->assign_solution(elite1->get_solution());
            elite1->generate_random();
        }
        for (i=0;i<population_size;i++)
        {
            fitnesses[i] = tabu_searchers[i]->get_fitness();
        }
        elite1_fitness = elite1->get_fitness();
        elite2_fitness = elite2->get_fitness();

        g++;
    }

    // PHASE 2 - ILS form of TabuCol with static tuned tabu tenure

    long cycles = 0;
    long cycles_max = ils_cycles_max;
    bool changed;

    if (! found)
    {

        printf("Current fitness: %lld", fitness);

        while (1)
        {
            if (cycles > cycles_max)
            {
                if (aux_output) printf("\n");
                break;
            }

            if (t > 0 && t % short_ils_length == 0)
            {
                if (aux_output)
                    if (changed)
                    {
                        printf("\nBest fitness: %lld (changed)\n", best_fitness);
                    }
                    else
                    {
                        printf("\nBest fitness: %lld\n", best_fitness);
                    }
                changed = false;
                tabu_searcher->assign_solution(best->get_solution());
                cycles++;
                previous_best_fitness = best_fitness;
            }

            fitness = tabu_searcher->get_fitness();

            if (0 == t % 1000)
            {
                if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
                {
                    putchar('\b');
                }
                if (aux_output) printf("%lld", fitness);
                old_fitness = fitness;
            }

            if (fitness >= G->m)
            {
                if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
                {
                    putchar('\b');
                }
                if (aux_output) printf("%lld\n", fitness);
                fprintf(f, "Found the following solution:\n");
                for (i=0;i<G->n;i++)
                {
                    fprintf(f, "%u", tabu_searcher->get_solution()[i]);
                    if (G->n-1 != i)
                    {
                        fprintf(f, ",");
                    }
                }
                fprintf(f, "\n");
                fprintf(f, "Iterations: %lld\n", t);
                fprintf(f, "Final number of conflicts: %lld\n", G->m-tabu_searcher->compute_fitness(tabu_searcher->get_solution()));
                found = true;
                break;
            }

            tabu_searcher->mutate_tabu();

            if (best_fitness <= tabu_searcher->get_fitness())
            {
                best_fitness = tabu_searcher->get_fitness();
                best->assign_solution(tabu_searcher->get_solution());
                changed = true;
            }

            t++;
        }
    }

    long improvement_cycles;

    if (! found)
    {
        for (improvement_cycles=0;improvement_cycles<improvement_cycles_max-1;improvement_cycles++)
        {
            if (found)
            {
                break;
            }

            // PHASE 3 - long-term H2Col with dynamic tabu tenure

            g = 0;
            for (i=0;i<population_size;i++)
            {
                if (improvement_cycles % 2 == 1)
                {
                    tabu_searchers[i]->init(G, colors, alpha, A, 0);
                }
                else
                {
                    tabu_searchers[i]->init(G, colors, 0, 0, B);
                }
                {
                    if (i == 0)
                    {
                        tabu_searchers[i]->assign_solution(best->get_solution());
                    }
                    else
                    {
                        tabu_searchers[i]->assign_solution(elite2->get_solution());
                    }
                    fitnesses[i] = tabu_searchers[i]->get_fitness();
                    for (j=0;j<ls_length;j++)
                    {
                        tabu_searchers[i]->mutate_tabu();
                    }
                    tabu_searchers[i]->assign_solution(tabu_searchers[i]->get_best_solution());
                }
                fitnesses[i] = tabu_searchers[i]->get_fitness();
            }

            if (improvement_cycles % 2 == 1)
            {
                offspring1->init(G, colors, alpha, A, 0);
            }
            else
            {
                offspring1->init(G, colors, 0, 0, B);
            }

            if (improvement_cycles % 2 == 1)
            {
                offspring2->init(G, colors, alpha, A, 0);
            }
            else
            {
                offspring2->init(G, colors, 0, 0, B);
            }

            if (improvement_cycles % 2 == 1)
            {
                elite1->init(G, colors, alpha, A, 0);
            }
            else
            {
                elite1->init(G, colors, 0, 0, B);
            }
            elite1->generate_random();
            elite1_fitness = elite1->get_fitness();

            if (improvement_cycles % 2 == 1)
            {
                elite2->init(G, colors, alpha, A, 0);
            }
            else
            {
                elite2->init(G, colors, 0, 0, B);
            }
            elite2->generate_random();
            elite2_fitness = elite2->get_fitness();

            stag = 0;
            while (1)
            {
                // fitnesses and their maximum
                max_fitness = 0;
                min_fitness = G->m;
                for (i=0;i<population_size;i++)
                {
                    if (max_fitness < fitnesses[i])
                    {
                        max_fitness_index = i;
                        max_fitness = fitnesses[i];
                    }
                    if (min_fitness > fitnesses[i])
                    {
                        min_fitness = fitnesses[i];
                    }
                }

                // auxiliary output
                if (aux_output) printf("Generation %lld: %lld, %lld; %lld. LS length: %ld\n", g, fitnesses[0], fitnesses[1], best_fitness, ls_length);

                // optimum handling
                if (max_fitness >= G->m)
                {
                    final_conflicts = G->m-tabu_searchers[max_fitness_index]->compute_fitness(tabu_searchers[max_fitness_index]->get_solution());
                    fprintf(f, "Found the following solution:\n");
                    for (i=0;i<G->n;i++)
                    {
                        fprintf(f, "%u", best->get_solution()[i]);
                        if (G->n-1 != i)
                        {
                            fprintf(f, ",");
                        }
                    }
                    fprintf(f, "\n");
                    fprintf(f, "Generations: %lld\n", g);
                    fprintf(f, "Iterations: %lld\n", t);
                    fprintf(f, "Final number of conflicts: %lld\n", G->m-best->compute_fitness(best->get_solution()));
                    found = true;

                    break;
                }

                // too close if at least 90% of the solution structure is the same
                similarity = partition_distance_heuristic(G, colors, tabu_searchers[0]->get_solution(), tabu_searchers[1]->get_solution());
                if (IMPOSSIBLE != similarity && similarity >= 99*G->n/100 || stag >= stag_max)
                {
                    final_conflicts = G->m-tabu_searchers[max_fitness_index]->compute_fitness(tabu_searchers[max_fitness_index]->get_solution());
                    break;
                }

                // selection of parents and crossover
                offspring1->create_by_greedy_crossover(tabu_searchers[0]->get_solution(), tabu_searchers[1]->get_solution());
                offspring2->create_by_greedy_crossover(tabu_searchers[1]->get_solution(), tabu_searchers[0]->get_solution());

                // tabu search on the offspring
                for (j=0;j<ls_length;j++)
                {
                    offspring1->mutate_tabu();
                    t++;
                }
                // tabu search on the offspring
                for (j=0;j<ls_length;j++)
                {
                    offspring2->mutate_tabu();
                    t++;
                }

                if (elite1_fitness <= offspring1->get_best_solution_fitness())
                {
                    elite1->assign_solution(offspring1->get_best_solution());
                    elite1_fitness = elite1->get_fitness();
                }
                if (elite1_fitness <= offspring2->get_best_solution_fitness())
                {
                    elite1->assign_solution(offspring2->get_best_solution());
                    elite1_fitness = elite1->get_fitness();
                }

                tabu_searchers[0]->assign_solution(offspring1->get_best_solution());
                tabu_searchers[1]->assign_solution(offspring2->get_best_solution());

                fitnesses[0] = tabu_searchers[0]->get_fitness();
                fitnesses[1] = tabu_searchers[1]->get_fitness();

                if (best_fitness <= elite1_fitness)
                {
                    best->assign_solution(elite1->get_solution());
                    best_fitness = elite1_fitness;
                    stag = 0;
                }
                else
                {
                    stag++;
                }

                if (g > 0 && g % iter_cycle == 0)
                {
                    tabu_searchers[0]->assign_solution(elite2->get_solution());
                    fitnesses[0] = tabu_searchers[0]->get_fitness();
                    elite2->assign_solution(elite1->get_solution());
                    elite1->generate_random();
                }
                for (i=0;i<population_size;i++)
                {
                    fitnesses[i] = tabu_searchers[i]->get_fitness();
                }
                elite1_fitness = elite1->get_fitness();
                elite2_fitness = elite2->get_fitness();

                g++;
            }

            // PHASE 4 - long-term ILS form of TabuCol with static tuned tabu tenure
            
            if (found)
            {
                break;
            }

            cycles = 0;
            // cycles_max = 100;

            printf("Current fitness: %lld", fitness);

            while (1)
            {
                if (cycles > cycles_max)
                {
                    if (aux_output) printf("\n");
                    break;
                }

                if (t > 0 && t % short_ils_length == 0)
                {
                    if (aux_output)
                        if (changed)
                        {
                            printf("\nBest fitness: %lld (changed)\n", best_fitness);
                        }
                        else
                        {
                            printf("\nBest fitness: %lld\n", best_fitness);
                        }
                    changed = false;
                    tabu_searcher->assign_solution(best->get_solution());
                    cycles++;
                    previous_best_fitness = best_fitness;
                }

                fitness = tabu_searcher->get_fitness();

                if (0 == t % 1000)
                {
                    if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
                    {
                        putchar('\b');
                    }
                    if (aux_output) printf("%lld", fitness);
                    old_fitness = fitness;
                }

                if (fitness >= G->m)
                {
                    if (aux_output) for (i=0;i<floor(log10(old_fitness))+1;i++)
                    {
                        putchar('\b');
                    }
                    if (aux_output) printf("%lld\n", fitness);
                    fprintf(f, "Found the following solution:\n");
                    for (i=0;i<G->n;i++)
                    {
                        fprintf(f, "%u", tabu_searcher->get_solution()[i]);
                        if (G->n-1 != i)
                        {
                            fprintf(f, ",");
                        }
                    }
                    fprintf(f, "\n");
                    fprintf(f, "Iterations: %lld\n", t);
                    fprintf(f, "Final number of conflicts: %lld\n", G->m-tabu_searcher->compute_fitness(tabu_searcher->get_solution()));
                    found = true;
                    break;
                }

                tabu_searcher->mutate_tabu();

                if (best_fitness <= tabu_searcher->get_fitness())
                {
                    best_fitness = tabu_searcher->get_fitness();
                    best->assign_solution(tabu_searcher->get_solution());
                    changed = true;
                }

                t++;
            }

        }
    }

    delete(tabu_searcher);

    for (i=0;i<population_size;i++)
    {
        delete(tabu_searchers[i]);
    }
    delete(offspring1);
    delete(offspring2);
    delete(elite1);
    delete(elite2);
    delete(best);

    *output_t = t;

    return (! found);
}
