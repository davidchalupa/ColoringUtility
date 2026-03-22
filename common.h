#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>

typedef unsigned int refer;

class vertex_list
{
public:
    bool occ;
    long vertex;
    vertex_list *next;
    vertex_list *prev;
    vertex_list *end;
};

class CommonSettings
{
public:
    static refer runs_greedy_clique;
    static refer runs_brelaz;
    static refer runs_rls_ig_clique;
    static refer runs_greedy_independent_set;
    static refer runs_rls_ig_independent_set;
    static refer runs_greedy_dominating_set;
    static refer runs_dfs;
    static refer runs_dfs_ls;
    static long long crossing_minimization_time_limit;
};

void QuickSort(refer *array, int left, int right);
int BinarySearch(refer *array, refer key, int left, int right);
unsigned long long power(unsigned long long base, unsigned long long exponent);

typedef struct
{
    // algorithm to use
    int algorithm;
    // generic parameters
    int agents_initial;
    int c;
    long t_max;
    // specific parameters
    int world_dim;
    int initial_food;
    int birth_period;
    int reap_ratio;
    int lifespan_ratio;
    int max_elite_list_size;
    int number_of_runs;
    int actual_run;
    int distance;
    int gamma;
    int crossovers;
    int tabu_list_size;
    int neighborhood;
    long long t_actual;
    long long max_fitness;
    unsigned long long init_time;
    // general data
    unsigned long long fitness_counter;
    unsigned long long mutation_counter;
    char elite_list_string[16384];
} algorithm_data;

algorithm_data *get_input_data();

#endif
