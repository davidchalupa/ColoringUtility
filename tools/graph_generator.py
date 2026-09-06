import networkx as nx

import numpy as np
from scipy.spatial import Delaunay

def generate_random_planar_graph(num_nodes):
    # random 2D points
    points = np.random.rand(num_nodes, 2)

    # compute Delaunay triangulation (always planar)
    tri = Delaunay(points)

    # create networkx graph from triangulation edges
    G = nx.Graph()
    for path in tri.simplices:
        G.add_edge(path[0], path[1])
        G.add_edge(path[1], path[2])
        G.add_edge(path[2], path[0])

    return G, points
