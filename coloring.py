import networkx as nx
import coloring_utility

#n = 100
#w = 4
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

n = 100000
w = 4
G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#G = nx.erdos_renyi_graph(100, 0.1, seed=42)


print(f"G: |V| = {G.number_of_nodes()}, |E| = {G.number_of_edges()}, d = {nx.density(G):.4f}")

try:
    colors = coloring_utility.process(G, time_limit=60)

    unique_colors = len(set(colors))

    print(f"Number of colors used: {unique_colors}")
    print("Colors:", colors)


except Exception as e:
    print(f"An error occurred in C++: {e}")
