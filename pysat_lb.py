import networkx as nx
from pysat.solvers import Glucose4

def solve_k_coloring(G, k):
    """
    Determines if a graph G is k-colorable using the Glucose4 SAT solver.
    Returns the coloring dictionary {node: color} if SAT, else None.
    """
    solver = Glucose4()
    nodes = list(G.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)

    # Helper to map (vertex_index, color_index) to a unique SAT variable ID
    # SAT IDs must be non-zero integers.
    def get_var(v_idx, c_idx):
        return v_idx * k + c_idx + 1

    # 1. Each vertex must have at least one color
    for i in range(n):
        solver.add_clause([get_var(i, c) for c in range(k)])

    # 2. Each vertex must have at most one color (prunes the search space)
    for i in range(n):
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                solver.add_clause([-get_var(i, c1), -get_var(i, c2)])

    # 3. Adjacent vertices cannot have the same color
    for u, v in G.edges():
        u_idx, v_idx = node_to_idx[u], node_to_idx[v]
        for c in range(k):
            solver.add_clause([-get_var(u_idx, c), -get_var(v_idx, c)])

    # Solve
    if solver.solve():
        model = solver.get_model()
        coloring = {}
        for i, node in enumerate(nodes):
            for c in range(k):
                if model[get_var(i, c) - 1] > 0:
                    coloring[node] = c
                    break
        return coloring
    else:
        return None

#n = 100
#w = 4
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#n = 100000
#w = 4
#G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

G = nx.erdos_renyi_graph(500, 0.1, seed=42)

k = 6
result = solve_k_coloring(G, k)
if result:
    print(f"Valid {k}-coloring found!")
else:
    print(f"No {k}-coloring exists.")
