#ifndef SAMPLING_H
#define SAMPLING_H

#include "tabu_base.h"
#include "graphs.h"

long long h2col(graph G, refer colors, int alpha, int A, int B, long long ls_length, long long *output_t, FILE *f);
long long newcol(graph G, refer colors, int alpha, int A, int B, long long ls_length, long short_ils_length, long ils_cycles_max, long improvement_cycles_max, long long *output_t, FILE *f);

#endif // SAMPLING_H
