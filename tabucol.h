#ifndef TABUCOL_H
#define TABUCOL_H

#include "tabu_base.h"
#include "graphs.h"

long long tabucol(graph G, refer colors, int alpha, int A, int B, long long ls_length, long long t_max, int stage, long long *output_t);

#endif // TABUCOL_H
