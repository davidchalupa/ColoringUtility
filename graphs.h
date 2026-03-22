#ifndef GRAPHS_H
#define GRAPHS_H
#include <stdio.h>
#include "common.h"

//#define MAX_VERTICES 5000050
#define MAX_VERTICES 50050
#define MAX_COLORS 1050
#define CLUSTERS_MAX 2000

typedef struct VERTEX
{
    refer edgecount;
    refer *sibl;
} vertex;

typedef struct GRAPH_DATA
{
    refer n;
    unsigned long m;
    double density;
    vertex V[MAX_VERTICES];
} graph_data;
typedef graph_data *graph;

class cluster
{
public:
    long genotype[MAX_VERTICES];
    long length;
    long long fitness;
    long defects;
    long parent;
};

int input_graph(FILE *source);
int output_graph(FILE *target);
void free_graph();
graph get_graph();
long get_problem();
bool are_adjacent(refer v, refer w);
void generate_graph_BA_model(unsigned long w, unsigned long n_max);
void generate_graph_UDG(unsigned long n_max, unsigned long range, unsigned long grid);
void generate_graph_WS_model(refer n_max, refer k_half, double beta);
void generate_graph_grid(refer rows, refer columns, double beta);
void generate_graph_complete_tree(refer branching_factor, refer depth);
void generate_complementary_graph(graph G);
void generate_graph_pruned_leaves(graph G);
void generate_largest_component(graph G);
void generate_shortcut_graph(graph G, refer k);
char *get_vertex_label(refer vertex);

#endif // GRAPHS_H
