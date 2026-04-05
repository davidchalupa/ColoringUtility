#ifndef TABUCOL_H
#define TABUCOL_H

#include "tabu_base.h"
#include "graphs.h"

#include <chrono>

long long tabucol(graph G, refer colors, int alpha, int A, int B, long long t_max, int stage, refer *result, std::chrono::high_resolution_clock::time_point &start_time, long long time_limit);

#endif // TABUCOL_H
