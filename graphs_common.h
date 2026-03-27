#ifndef GRAPHS_COMMON_H
#define GRAPHS_COMMON_H

#include "common.h"

#define MAX_VERTICES 5000050

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


#endif // GRAPHS_COMMON_H
