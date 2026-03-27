import networkx as nx
import coloring_utility

G = nx.Graph()
G.add_edges_from([
    ("A", "B"),
    ("B", "C"),
    ("C", "A"),
    ("C", "D")
])

print("Passing graph to C++...\n")

coloring_utility.process(G)

print("\nSuccessfully returned to Python!")
