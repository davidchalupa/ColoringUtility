#ifndef CLI_H
#define CLI_H
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "graphs.h"
#include "random_generator.h"
#include "time.h"
#include "statistics.h"

#define INSTANCE_LONG_CYCLE_FIXED 1
#define INSTANCE_LONG_CYCLE 2
#define INSTANCE_HAMILTONIAN_CYCLE 3

class cli
{
private:
    char filename[PATH_MAX];
    char filename_output[PATH_MAX];
    FILE *source;
    graph G;
    // statistics
    unsigned long max_deg,min_deg,triangles;
    unsigned long long current_time,avg_time,avg_iter;
    double avg_deg,stdev_deg;
    refer degree_distrib[MAX_VERTICES];
    // instance parameters
    refer k;
    unsigned long w,n_max,range,grid;
    // other
    refer i;
    int tabu_tenure;
    unsigned long long t_max;
    unsigned long long ls_length;
    long ils_cycles_max;
    long improvement_cycles_max;
    refer runs;
    // problem-specific
    int choose_instance();
    int generate_instance();
    int generate_instance_udg();
    int compute_statistics();
    // management methods
    void usage();
    void process_params();
    void sleep(int milisec);
    void compute();
public:
    cli();
    int start_cli(int argc, char **argv);
};

#endif // CLI_H
