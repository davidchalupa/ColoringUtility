import networkx as nx
import coloring_utility

n = 100000
w = 4
G = nx.barabasi_albert_graph(n=n, m=w)

print(f"G: |V| = {G.number_of_nodes()}, |E| = {G.number_of_edges()}, d = {nx.density(G):.4f}")

try:
    coloring_utility.process(G, time_limit=600)
except Exception as e:
    print(f"An error occurred in C++: {e}")
