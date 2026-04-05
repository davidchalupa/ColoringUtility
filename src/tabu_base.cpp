#include <math.h>

#include "tabu_base.h"

#include "algorithm.h"
#include "algorithm_brelaz.h"

tabu_base::tabu_base()
{
}

tabu_base::~tabu_base()
{
}

void tabu_base::deinit()
{
    delete[](solution);
    delete[](best_solution);
    delete[](tabu_list);

    for (refer v = 0; v < colors; v++)
    {
        delete[](tabu_matrix[v]);
    }
    delete[](tabu_matrix);
    for (refer v = 0; v < colors; v++)
    {
        delete[](neighbor_color_matrix[v]);
    }
    delete[](neighbor_color_matrix);

    delete[](conflicting_vertices_list);
    delete[](positioning_list);
    delete[](collisions);
    delete[](best_vertices);
    delete[](best_colors);
    delete[](clusters);
}

void tabu_base::init(graph G, refer colors, int alpha, int A, int B)
{
    refer v,c;

    solution = new refer[G->n];
    best_solution = new refer[G->n];
    tabu_list = new list[MAX_TABU_LIST_SIZE];

    tabu_matrix = new bool*[G->n];
    for (refer v = 0; v < G->n; v++)
    {
        tabu_matrix[v] = new bool[colors+1];
    }
    neighbor_color_matrix = new short*[G->n];
    for (refer v = 0; v < G->n; v++)
    {
        neighbor_color_matrix[v] = new short[colors+1];
    }

    conflicting_vertices_list = new long[G->n];
    positioning_list = new long[G->n];
    collisions = new long[G->n];
    best_vertices = new refer[G->n*(colors+1)];
    best_colors = new refer[G->n*(colors+1)];
    clusters = new cluster[CLUSTERS_MAX];

    this->colors = colors;
    this->G = G;
    this->alpha = alpha;
    this->A = A;
    this->B = B;

    tabu_list_size = 0;
    max_tabu_list_size = MAX_TABU_LIST_SIZE;

    fitness = 0;
    fitness_counter = 0;

    for (v=0;v<G->n;v++)
    {
        for (c=1;c<=colors;c++)
        {
            tabu_matrix[v][c] = false;
            neighbor_color_matrix[v][c] = 0;
        }
    }

    for (v=0;v<G->n;v++)
    {
        collisions[v] = 0;
        positioning_list[v] = 0;
        conflicting_vertices_list[v] = 0;
        best_vertices[v] = 0;
        best_colors[v] = 0;
        solution[v] = 0;
        best_solution[v] = 0;
    }

    tabu_list_begin = NULL;
    previous_vertex = 0;
    previous_color = 0;
    previous_fitness = 0;
    aspiration_fitness = 0;
    best_moves_count = 0;
    conflicts = 0;
    fitness = 0;
    best_state_fitness = 0;

}

long long tabu_base::compute_fitness(refer *solution)
{
    refer i,v,c;
    long long cnt = G->m;

    // deleting the neighborhood color matrix
    for (v=0;v<G->n;v++)
    {
        for (c=1;c<=colors;c++)
        {
            neighbor_color_matrix[v][c] = 0;
        }
    }

    // deleting the collisions vector
    for (v=0;v<G->n;v++)
    {
        collisions[v] = 0;
    }

    for (v=0;v<G->n;v++)
    {
        for (i=0;i<G->V[v].edgecount;i++)
        {
            if (G->V[v].sibl[i] > v)
            {
                // setting the matrix
                neighbor_color_matrix[v][solution[G->V[v].sibl[i]]]++;
                neighbor_color_matrix[G->V[v].sibl[i]][solution[v]]++;
                // setting the collisions and fitness
                if (solution[G->V[v].sibl[i]] == solution[v])
                {
                    cnt--;
                    collisions[v]++;
                    collisions[G->V[v].sibl[i]]++;
                }
            }
        }
    }

    conflicts = 0;
    for (v=0;v<G->n;v++)
    {
        if (collisions[v] != 0)
        {
            conflicting_vertices_list[conflicts] = v;
            positioning_list[v] = conflicts;
            conflicts++;
        }
    }

    fitness_counter++;

    return cnt;
}

void tabu_base::generate_random()
{
    refer i;

    for (i=0;i<G->n;i++)
    {
        solution[i] = generator.random(1,colors);
    }
    /*solution[0] = 1;
    for (i=0;i<(G->n-1)/5;i++)
    {
        if (i % 3 == 0)
        {
            solution[5*i+1] = 2;
            solution[5*i+2] = 3;
            solution[5*i+3] = 1;
            solution[5*i+4] = 1;
            solution[5*i+5] = 1;
        }
        else if (i % 3 == 1)
        {
            solution[5*i+1] = 1;
            solution[5*i+2] = 2;
            solution[5*i+3] = 3;
            solution[5*i+4] = 3;
            solution[5*i+5] = 3;
        }
        else
        {
            solution[5*i+1] = 1;
            solution[5*i+2] = 3;
            solution[5*i+3] = 2;
            solution[5*i+4] = 2;
            solution[5*i+5] = 2;
        }
    }*/

    for (i=0;i<G->n;i++)
    {
        best_solution[i] = solution[i];
    }

    fitness = compute_fitness(solution);

    aspiration_fitness = fitness;
    // ToDo: maybe change
    //best_state_fitness = fitness;
    best_state_fitness = 0;
}

void tabu_base::generate_Brelaz()
{
    refer i;

    algorithm *algorithm_brelaz_instance = new algorithm_brelaz();

    algorithm_brelaz_instance->brelaz_with_heap(G, solution);

    for (i=0;i<G->n;i++)
    {
        if (solution[i] > colors)
        {
            solution[i] = generator.random(1,colors);
        }
    }

    delete(algorithm_brelaz_instance);

    for (i=0;i<G->n;i++)
    {
        best_solution[i] = solution[i];
    }

    fitness = compute_fitness(solution);

    aspiration_fitness = fitness;
    // ToDo: maybe change
    //best_state_fitness = fitness;
    best_state_fitness = 0;
}

/////////////////
// TABU SEARCH //
/////////////////

bool tabu_base::is_tabu(long vertex, long color)
{
    return tabu_matrix[vertex][color];
}

void tabu_base::decrease_tabu()
{
    long j;

    // we decrease the tabu tenures of existing tabu moves
    for (j=0;j<tabu_list_size;j++)
    {
        tabu_list[j].tabu_tenure--;
    }
    j = 0;
    while (j < tabu_list_size)
    {
        // if we reached zero, we erase the move
        if (tabu_list[j].tabu_tenure == 0)
        {
            tabu_matrix[tabu_list[j].vertex][tabu_list[j].color] = false;
            tabu_list[j].vertex = tabu_list[tabu_list_size-1].vertex;
            tabu_list[j].color = tabu_list[tabu_list_size-1].color;
            tabu_list[j].tabu_tenure = tabu_list[tabu_list_size-1].tabu_tenure;
            tabu_list_size--;
        }
        else
        {
            tabu_matrix[tabu_list[j].vertex][tabu_list[j].color] = true;
            j++;
        }
    }
}

void tabu_base::insert_tabu(long vertex, long color, long conflicts)
{
    // we insert the move
    if (tabu_list_size < max_tabu_list_size /*&& (alpha > 0 || A > 0)*/)
    {
        // dynamic tabu tenure with a static component
        tabu_list[tabu_list_size].tabu_tenure = alpha * conflicts / 10 + generator.random(0, A-1) + B;
        if (tabu_list[tabu_list_size].tabu_tenure <= 0)
        {
            return;
        }

        tabu_matrix[vertex][color] = true;
        tabu_list[tabu_list_size].vertex = vertex;
        tabu_list[tabu_list_size].color = color;
        if (tabu_list[tabu_list_size].tabu_tenure < 1)
        {
            tabu_list[tabu_list_size].tabu_tenure = 1;
        }
        tabu_list_size++;
    }
}

long tabu_base::mutate_state_stochastic(long vertex, long c, long initial_fitness, long old_color)
{
    return initial_fitness + neighbor_color_matrix[vertex][old_color] - neighbor_color_matrix[vertex][c];
}

void tabu_base::update_collisions(long v, refer old_color)
{
    refer i,new_color;
    long old_collisions,old_collisions_center;
    new_color = solution[v];

    old_collisions_center = collisions[v];
    for (i=0;i<G->V[v].edgecount;i++)
    {
        old_collisions = collisions[G->V[v].sibl[i]];
        // if the edge was changed from monochromatic to correct
        if (solution[G->V[v].sibl[i]] == old_color && solution[G->V[v].sibl[i]] != new_color)
        {
            collisions[v]--;
            collisions[G->V[v].sibl[i]]--;
        }
        // if the edge was changed from correct to monochromatic
        if (solution[G->V[v].sibl[i]] == new_color && solution[G->V[v].sibl[i]] != old_color)
        {
            collisions[v]++;
            collisions[G->V[v].sibl[i]]++;
        }
        // if we resolved a local conflict
        if (collisions[G->V[v].sibl[i]] == 0 && old_collisions != 0)
        {
            conflicting_vertices_list[positioning_list[G->V[v].sibl[i]]] = conflicting_vertices_list[conflicts-1];
            positioning_list[conflicting_vertices_list[conflicts-1]] = positioning_list[G->V[v].sibl[i]];
            positioning_list[G->V[v].sibl[i]] = -1;
            conflicts--;
        }
        // if we created a local conflict
        if (collisions[G->V[v].sibl[i]] != 0 && old_collisions == 0)
        {
            conflicting_vertices_list[conflicts] = G->V[v].sibl[i];
            positioning_list[G->V[v].sibl[i]] = conflicts;
            conflicts++;
        }
        // we change the neighborhood color matrix
        neighbor_color_matrix[G->V[v].sibl[i]][old_color]--;
        neighbor_color_matrix[G->V[v].sibl[i]][new_color]++;
    }
    // if we resolved a conflict on the center vertex
    if (collisions[v] == 0 && old_collisions_center != 0)
    {
        conflicting_vertices_list[positioning_list[v]] = conflicting_vertices_list[conflicts-1];
        positioning_list[conflicting_vertices_list[conflicts-1]] = positioning_list[v];
        positioning_list[v] = -1;
        conflicts--;
    }
    // if we created a conflict on the center vertex
    if (collisions[v] != 0 && old_collisions_center == 0)
    {
        conflicting_vertices_list[conflicts] = v;
        positioning_list[v] = conflicts;
        conflicts++;
    }

}

void tabu_base::mutate_tabu()
{
    refer i;
    long j,r;
    refer c;
    refer actual_color;
    long new_fitness;
    long best_fitness;
    long vertex,color;
    long index;

    // we construct the neighborhood
    best_fitness = 0;
    best_moves_count = 0;
    for (r=0;r<conflicts;r++)
    {
        j = conflicting_vertices_list[r];
        actual_color = solution[j];
        for (c=1;c<=colors;c++)
        {
            if (c == actual_color)
            {
                continue;
            }

            // we calculate the potential new fitness
            new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

            // if the move is good, we keep the information
            if (best_fitness <= new_fitness && (! is_tabu(j,c) || aspiration_fitness < new_fitness))
            {
                if (best_fitness == new_fitness)
                {
                    best_vertices[best_moves_count] = j;
                    best_colors[best_moves_count] = c;
                    best_moves_count++;
                }
                else
                {
                    best_moves_count = 0;
                    best_vertices[best_moves_count] = j;
                    best_colors[best_moves_count] = c;
                    best_moves_count++;
                    best_fitness = new_fitness;
                }
            }
        }
    }
    fitness_counter+=conflicts*(colors-1);
    //iteration_counter++;

    if (best_fitness > fitness)
    {
        locopt = false;
    }
    else
    {
        locopt = true;
    }

    decrease_tabu();

    if (0 == best_moves_count)
    {
        return;
    }

    index = generator.random(0,best_moves_count-1);
    vertex = best_vertices[index];
    color = best_colors[index];

    bool new_best_solution = false;
    if (aspiration_fitness <= best_fitness)
    {
        new_best_solution = true;
        aspiration_fitness = best_fitness;
    }

    previous_color = solution[vertex];
    insert_tabu(vertex,previous_color,conflicts);

    if (fitness == best_fitness)
    {
        t_stag++;
    }
    else
    {
        t_stag = 0;
    }
    fitness = best_fitness;
    solution[vertex] = color;
    if (new_best_solution)
    {
        for (i=0;i<G->n;i++)
        {
            best_solution[i] = solution[i];
        }
        best_state_fitness = best_fitness;
    }
    update_collisions(vertex,previous_color);
}

void tabu_base::mutate_tabu_rls_b()
{
    refer i;
    long j,r;
    refer c;
    refer actual_color;
    long new_fitness;
    long best_fitness;
    long vertex,color;
    refer index;

    // we construct the neighborhood
    best_fitness = 0;
    best_moves_count = 0;
    for (r=0;r<conflicts;r++)
    {
        j = conflicting_vertices_list[r];
        actual_color = solution[j];
        for (c=1;c<=colors;c++)
        {
            if (c == actual_color)
            {
                continue;
            }

            // we calculate the potential new fitness
            new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

            // if the move is good, we keep the information
            if (fitness < new_fitness && (! is_tabu(j,c) || aspiration_fitness < new_fitness))
            {
                best_vertices[best_moves_count] = j;
                best_colors[best_moves_count] = c;
                best_moves_count++;
            }
        }
    }
    fitness_counter+=conflicts*(colors-1);
    //iteration_counter++;

    if (best_moves_count)
    {
        index = generator.random(0,best_moves_count-1);
        vertex = best_vertices[index];
        color = best_colors[index];
        actual_color = solution[vertex];

        new_fitness = mutate_state_stochastic(vertex,color,fitness,actual_color);

        locopt = false;
    }
    else
    {
        do
        {
            r = generator.random(0,conflicts-1);
            j = conflicting_vertices_list[r];
            actual_color = solution[j];
            c = generator.random(1,colors);
        }
        while (c == actual_color);

        vertex = j;
        color = c;
        new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

        locopt = true;
    }

    best_fitness = new_fitness;

    decrease_tabu();
    bool new_best_solution = false;
    if (aspiration_fitness <= best_fitness)
    {
        new_best_solution = true;
        aspiration_fitness = best_fitness;
    }

    previous_color = solution[vertex];
    insert_tabu(vertex,previous_color,conflicts);

    if (fitness == new_fitness)
    {
        t_stag++;
    }
    else
    {
        t_stag = 0;
    }
    fitness = new_fitness;
    solution[vertex] = color;
    if (new_best_solution)
    {
        for (i=0;i<G->n;i++)
        {
            best_solution[i] = solution[i];
        }
        best_state_fitness = best_fitness;
    }
    update_collisions(vertex,previous_color);
}

void tabu_base::mutate_tabu_rls()
{
    refer i;
    long j,r;
    refer c;
    refer actual_color;
    long new_fitness;
    long best_fitness;
    long vertex,color;
    long index;

    best_moves_count = 0;

    locopt = false;

    for (i=0;i<conflicts*(colors-1);i++)
    {
        do
        {
            r = generator.random(0,conflicts-1);
            j = conflicting_vertices_list[r];
            actual_color = solution[j];
            c = generator.random(1,colors);
        }
        while (c == actual_color);

        new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

        // if the move is good, we keep the information
        if (fitness < new_fitness && (! is_tabu(j,c) || aspiration_fitness < new_fitness))
        {
            vertex = j;
            color = c;
            best_moves_count = 1;
            break;
        }
        fitness_counter++;
    }
    //iteration_counter++;

    if (! best_moves_count)
    {
        do
        {
            r = generator.random(0,conflicts-1);
            j = conflicting_vertices_list[r];
            actual_color = solution[j];
            c = generator.random(1,colors);
        }
        while (c == actual_color);

        vertex = j;
        color = c;
        new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

        locopt = true;
    }

    best_fitness = new_fitness;

    decrease_tabu();
    bool new_best_solution = false;
    if (aspiration_fitness <= best_fitness)
    {
        new_best_solution = true;
        aspiration_fitness = best_fitness;
    }

    previous_color = solution[vertex];
    insert_tabu(vertex,previous_color,conflicts);

    if (fitness == new_fitness)
    {
        t_stag++;
    }
    else
    {
        t_stag = 0;
    }
    fitness = new_fitness;
    solution[vertex] = color;
    if (new_best_solution)
    {
        for (i=0;i<G->n;i++)
        {
            best_solution[i] = solution[i];
        }
        best_state_fitness = best_fitness;
    }
    update_collisions(vertex,previous_color);
}

/*void tabu_base::mutate_tabu_rls()
{
    refer i;
    long j,r;
    refer c;
    refer actual_color;
    long new_fitness;
    long best_fitness;
    long vertex,color;
    long index;

    best_moves_count = 0;

    locopt = false;

    for (i=0;i<conflicts*(colors-1);i++)
    {
        do
        {
            r = generator.random(0,conflicts-1);
            j = conflicting_vertices_list[r];
            actual_color = solution[j];
            c = generator.random(1,colors);
        }
        while (c == actual_color);

        new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

        // if the move is good, we keep the information
        if (fitness < new_fitness && (! is_tabu(j,c) || aspiration_fitness < new_fitness))
        {
            vertex = j;
            color = c;
            best_moves_count = 1;
            break;
        }
        fitness_counter++;
    }
    //iteration_counter++;

    if (! best_moves_count)
    {
        do
        {
            r = generator.random(0,conflicts-1);
            j = conflicting_vertices_list[r];
            actual_color = solution[j];
            c = generator.random(1,colors);
        }
        while (c == actual_color);

        locopt = true;
    }

    best_fitness = fitness;

    decrease_tabu();
    bool new_best_solution = false;
    if (aspiration_fitness <= best_fitness)
    {
        new_best_solution = true;
        aspiration_fitness = best_fitness;
    }

    previous_color = solution[vertex];
    insert_tabu(vertex,previous_color,conflicts);

    if (fitness == new_fitness)
    {
        t_stag++;
    }
    else
    {
        t_stag = 0;
    }
    fitness = new_fitness;
    solution[vertex] = color;
    if (new_best_solution)
    {
        for (i=0;i<G->n;i++)
        {
            best_solution[i] = solution[i];
        }
        best_state_fitness = best_fitness;
    }
    update_collisions(vertex,previous_color);
}*/

bool tabu_base::is_probably_local_optimum()
{
    return locopt;
}

void tabu_base::mutate_random()
{
    refer i,c;
    long j,r;
    long new_fitness;
    long vertex,color;
    long best_fitness;
    refer actual_color;

    // we construct the neighborhood
    r = generator.random(0, conflicts-1);
    j = conflicting_vertices_list[r];
    actual_color = solution[j];
    c = generator.random(1, colors);
    // we calculate the potential new fitness
    new_fitness = mutate_state_stochastic(j,c,fitness,actual_color);

    fitness_counter+=1;
    //iteration_counter++;

    decrease_tabu();

    vertex = j;
    color = c;
    best_fitness = new_fitness;

    bool new_best_solution = false;
    if (aspiration_fitness <= best_fitness)
    {
        new_best_solution = true;
        aspiration_fitness = best_fitness;
    }

    previous_color = solution[vertex];
    insert_tabu(vertex,previous_color,conflicts);

    fitness = best_fitness;
    solution[vertex] = color;
    if (new_best_solution)
    {
        for (i=0;i<G->n;i++)
        {
            best_solution[i] = solution[i];
        }
        best_state_fitness = best_fitness;
    }
    update_collisions(vertex,previous_color);
}

refer *tabu_base::get_solution()
{
    return solution;
}

refer *tabu_base::get_best_solution()
{
    return best_solution;
}

long long tabu_base::get_best_solution_fitness()
{
    return best_state_fitness;
}

long long tabu_base::get_fitness()
{
    return fitness;
}

void tabu_base::set_alpha_A(int alpha, int A)
{
    this->alpha = alpha;
    this->A = A;
}

void tabu_base::assign_solution(refer *solution)
{
    refer i;
    refer auxiliary[MAX_VERTICES_TABU_BASE];

    for (i=0;i<G->n;i++)
    {
        auxiliary[i] = solution[i];
    }
    init(G, colors, alpha, A, B);

    for (i=0;i<G->n;i++)
    {
        this->solution[i] = auxiliary[i];
        this->best_solution[i] = this->solution[i];
    }

    fitness = compute_fitness(this->solution);

    aspiration_fitness = fitness;
    // ToDo: maybe change
    //best_state_fitness = fitness;
    best_state_fitness = 0;
}

void tabu_base::assign_solution_keep_aspiration(refer *solution)
{
    refer i;
    refer auxiliary[MAX_VERTICES_TABU_BASE];

    for (i=0;i<G->n;i++)
    {
        auxiliary[i] = solution[i];
    }

    for (i=0;i<G->n;i++)
    {
        this->solution[i] = auxiliary[i];
    }

    fitness = compute_fitness(this->solution);
}

///////////////
// CROSSOVER //
///////////////

long tabu_base::extract_cluster(refer *parent, long *cluster, refer r)
{
    // we seek for all nodes with chosen color
    refer i;
    long j;

    j = 0;
    for (i=0;i<G->n;i++)
    {
        // if we find the matching color, we insert it
        if (parent[i] == r)
        {
            cluster[j] = i;
            j++;
        }
    }

    return j;
}

long tabu_base::receive_cluster(long *cluster, long length, long color)
{
    long i;

    // we assign each color
    for (i=0;i<length;i++)
    {
        solution[cluster[i]] = color;
    }

    return 0;
}

long tabu_base::change_fitness_after_mutation(long v, long initial_fitness, long old_color)
{
    long new_fitness,new_color;

    new_fitness = initial_fitness;
    new_color = solution[v];

    new_fitness += neighbor_color_matrix[v][old_color];
    new_fitness -= neighbor_color_matrix[v][new_color];

    return new_fitness;
}

long tabu_base::mutate_state_augment(long vertex, long initial_fitness, bool skip_former_color)
{
    refer c;
    long best_color;
    long best_collisions;
    refer old_color;
    long fitness;

    old_color = solution[vertex];

    // we find the minimal collisioning color
    best_color = 1;
    best_collisions = G->n;
    for (c=1;c<=colors;c++)
    {
        if (skip_former_color && c == old_color)
        {
            continue;
        }
        if (best_collisions > neighbor_color_matrix[vertex][c])
        {
            best_color = c;
            best_collisions = neighbor_color_matrix[vertex][best_color];
        }
    }

    solution[vertex] = best_color;
    fitness = change_fitness_after_mutation(vertex,initial_fitness,old_color);
    update_collisions(vertex,old_color);

    return fitness;
}

long tabu_base::get_conflicting_vertices_count()
{
    return conflicts;
}

void tabu_base::create_by_greedy_crossover(refer *parent1, refer *parent2)
{
    long i,j,k,l,t,c,clusters_count;
    long max_cluster_fitness_index = 0;

    init(G, colors, alpha, A, B);

    clusters_count = 0;
    for (c=1;c<=colors;c++)
    {
        clusters[clusters_count].length = extract_cluster(parent1,clusters[clusters_count].genotype,c);
        clusters[clusters_count].fitness = clusters[clusters_count].length;
        clusters[clusters_count].parent = 0;
        clusters_count++;
        clusters[clusters_count].length = extract_cluster(parent2,clusters[clusters_count].genotype,c);
        clusters[clusters_count].fitness = clusters[clusters_count].length;
        clusters[clusters_count].parent = 1;
        clusters_count++;
    }

    long long max_cluster_fitness = 0;

    // no used clusters
    bool used_clusters[CLUSTERS_MAX];
    for (i=0;i<clusters_count;i++)
    {
        used_clusters[i] = false;
    }

    int current_parent = 0;

    t = 0;
    while (1)
    {
        max_cluster_fitness = 0;
        for (k=0;k<clusters_count;k++)
        {
            if (! used_clusters[k] && clusters[k].parent == current_parent && max_cluster_fitness < clusters[k].fitness)
            {
                max_cluster_fitness = clusters[k].fitness;
                max_cluster_fitness_index = k;
            }
        }

        used_clusters[max_cluster_fitness_index] = true;
        for (j=0;j<clusters[max_cluster_fitness_index].length;j++)
        {
            // erasing the vertex from all other clusters
            for (k=0;k<clusters_count;k++)
            {
                if (k == max_cluster_fitness_index)
                {
                    continue;
                }
                for (l=0;l<clusters[k].length;l++)
                {
                    if (clusters[k].genotype[l] == clusters[max_cluster_fitness_index].genotype[j])
                    {
                        clusters[k].genotype[l] = clusters[k].genotype[clusters[k].length-1];
                        clusters[k].length--;
                        clusters[k].fitness = clusters[k].length;
                    }
                }
            }
        }
        receive_cluster(clusters[max_cluster_fitness_index].genotype,clusters[max_cluster_fitness_index].length,t+1);
        t++;

        if (t >= colors)
        {
            break;
        }
        current_parent++;
        current_parent %= 2;
    }

    long q = 0;
    for (i=0;i<G->n;i++)
    {
        if (0 == solution[i])
        {
            solution[i] = generator.random(1, colors);
            q++;
        }
    }

    fitness = compute_fitness(solution);

    for (i=0;i<G->n;i++)
    {
        best_solution[i] = solution[i];
    }

    aspiration_fitness = fitness;
    best_state_fitness = fitness;
}

