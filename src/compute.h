#ifndef COMPUTE_H
#define COMPUTE_H

#include "common.h"

#include "graphs_common.h"

void compute(graph G, refer *coloring, refer &best_lower_bound, long long time_limit, bool apply_sat_lower_bound);

#endif // COMPUTE_H
