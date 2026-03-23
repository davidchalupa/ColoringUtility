#ifndef TABU_BASE_H
#define TABU_BASE_H
#include "graphs.h"
#include "random_generator.h"

#define MAX_TABU_LIST_SIZE 5000

#define MAX_VERTICES_TABU_BASE 50050

#define CLUSTERS_MAX 2000

class cluster
{
public:
    long genotype[MAX_VERTICES];
    long length;
    long long fitness;
    long defects;
    long parent;
};

class tabu_base
{
private:
    class list
    {
    public:
        long vertex;
        long color;
        long tabu_tenure;
        list *next;
    };

    graph G;
    refer colors;
    long long fitness;

    list *tabu_list_begin;
    long tabu_list_size;
    long max_tabu_list_size;

//    refer solution[MAX_VERTICES_TABU_BASE];
//    refer best_solution[MAX_VERTICES_TABU_BASE];
//    list tabu_list[MAX_TABU_LIST_SIZE];
//    bool tabu_matrix[MAX_VERTICES_TABU_BASE][MAX_COLORS];
//    short neighbor_color_matrix[MAX_VERTICES_TABU_BASE][MAX_COLORS];
//    long conflicting_vertices_list[MAX_VERTICES_TABU_BASE];
//    long positioning_list[MAX_VERTICES_TABU_BASE];
//    long collisions[MAX_VERTICES_TABU_BASE];
//    short best_vertices[MAX_VERTICES_TABU_BASE*MAX_COLORS];
//    short best_colors[MAX_VERTICES_TABU_BASE*MAX_COLORS];
//    cluster clusters[CLUSTERS_MAX];

    refer *solution;
    refer *best_solution;
    list *tabu_list;
    bool **tabu_matrix;
    short **neighbor_color_matrix;
    long *conflicting_vertices_list;
    long *positioning_list;
    long *collisions;
    refer *best_vertices;
    refer *best_colors;
    cluster *clusters;

    long aspiration_fitness;
    long long best_state_fitness;
    long conflicts;
    long previous_vertex;
    long previous_color;
    long previous_fitness;
    long best_moves_count;

    bool locopt;

    refer alpha;
    refer A;
    refer B;
    unsigned long long t_stag;

    long long fitness_counter;

    random_generator generator;

    long mutate_state_stochastic(long vertex, long c, long initial_fitness, long old_color);
    void update_collisions(long v, refer old_color);
    bool is_tabu(long vertex, long color);
    void decrease_tabu();
    void insert_tabu(long vertex, long color, long tenure);
    long extract_cluster(refer *parent, long *cluster, refer r);
    long receive_cluster(long cluster[], long length, long color);
    long mutate_state_augment(long vertex, long initial_fitness, bool skip_former_color);
    long change_fitness_after_mutation(long v, long initial_fitness, long old_color);
public:
    tabu_base();
    ~tabu_base();
    void init(graph G, refer colors, int alpha, int A, int B);
    void deinit();
    void generate_random();
    void generate_Brelaz();
    long long compute_fitness(refer *solution);
    long long get_fitness();
    refer *get_solution();
    refer *get_best_solution();
    long long get_best_solution_fitness();
    void mutate_tabu();
    void mutate_tabu_rls();
    void mutate_tabu_rls_b();
    void mutate_random();
    void set_alpha_A(int alpha, int A);
    void assign_solution(refer *solution);
    void assign_solution_keep_aspiration(refer *solution);
    long get_conflicting_vertices_count();
    void create_by_greedy_crossover(refer *parent1, refer *parent2);
    bool is_probably_local_optimum();
};

#endif // TABU_BASE_H
