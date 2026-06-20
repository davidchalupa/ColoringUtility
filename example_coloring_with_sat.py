import networkx as nx
import coloring_utility

from tools import coloring_sat

#n = 100
#w = 3
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#n = 100
#w = 4
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

n = 100
w = 5
G = nx.barabasi_albert_graph(n=n, m=w, seed=442)

#n = 100000
#w = 4
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#G = nx.erdos_renyi_graph(100, 0.1, seed=42)


print(f"G: |V| = {G.number_of_nodes()}, |E| = {G.number_of_edges()}, d = {nx.density(G):.4f}")

try:
    colors, lower_bound = coloring_utility.process(G, time_limit=60)

    unique_colors = len(set(colors))

    print(f"Number of colors used: {unique_colors}")
    print(f"Best lower bound found: {lower_bound}")

    if lower_bound < unique_colors:
        print("Optimal solution not yet found. Attempting to close the gap with SAT...")
        k = lower_bound
        while k < unique_colors:
            result = coloring_sat.solve_k_coloring(G, k)
            if result:
                print(f"Valid {k}-coloring found!")
                # ToDo: override the result from the C++ solver if SAT was better
                break
            else:
                print(f"No {k}-coloring exists.")
                k += 1
                lower_bound = k
        print(f"Improved lower bound: {lower_bound}")

        if lower_bound == unique_colors:
            print(f"Found an optimal coloring with {unique_colors} colors!")

    print("Colors:", colors)

except Exception as e:
    print(f"An error occurred in C++: {e}")

