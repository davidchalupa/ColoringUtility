import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import time

def compute_fractional_chromatic_bound(G, max_iter=150, time_limit=300):
    """
    Computes a tight lower bound for the fractional chromatic number
    using Column Generation.
    """
    n = G.number_of_nodes()
    if n == 0:
        return 0.0

    nodes = list(G.nodes())
    node_idx = {v: i for i, v in enumerate(nodes)}

    # ---------------------------------------------------------
    # 1. Setup the Pricing Problem Constraints (Sparse Matrix)
    # ---------------------------------------------------------
    # Constraint: y_u + y_v <= 1 for all edges (u,v)
    edges = list(G.edges())
    num_edges = len(edges)

    row = []
    col = []
    data = []
    for i, (u, v) in enumerate(edges):
        row.extend([i, i])
        col.extend([node_idx[u], node_idx[v]])
        data.extend([1, 1])

    A_pricing = sp.csr_matrix((data, (row, col)), shape=(num_edges, n))
    b_pricing = np.ones(num_edges)

    pricing_constraints = LinearConstraint(A_pricing, ub=b_pricing)
    bounds = Bounds(0, 1)
    integrality = np.ones(n) # All variables must be integers (0 or 1)

    # ---------------------------------------------------------
    # 2. Initialize Master Problem
    # ---------------------------------------------------------
    # Start with n trivial independent sets (each containing exactly 1 vertex)
    # A_master is a (nodes x independent_sets) matrix
    A_master = np.eye(n)

    best_lower_bound = 0.0
    start_time = time.time()

    print(f"Starting Column Generation for Graph ({n} nodes, {num_edges} edges)")
    print("-" * 65)
    print(f"{'Iter':<5} | {'Z_RMP (Upper)':<15} | {'Max Weight (W)':<15} | {'Lower Bound':<15}")
    print("-" * 65)

    for iteration in range(max_iter):
        if time.time() - start_time > time_limit:
            print("\nTime limit reached. Halting early.")
            break

        # ---------------------------------------------------------
        # 3. Solve the Dual of the Restricted Master Problem
        # ---------------------------------------------------------
        # Dual Variables: pi_v for v in V.
        # Maximize sum(pi_v) subject to A_master^T * pi <= 1
        # SciPy minimizes, so we minimize -sum(pi_v)
        c_dual = -np.ones(n)
        A_ub_dual = A_master.T
        b_ub_dual = np.ones(A_master.shape[1])

        res_dual = linprog(c_dual, A_ub=A_ub_dual, b_ub=b_ub_dual,
                           bounds=(0, None), method='highs')

        if not res_dual.success:
            raise RuntimeError("Master Dual LP failed to converge.")

        Z_rmp = -res_dual.fun # Current fractional coloring value (Upper Bound)
        pi = res_dual.x       # Dual values (vertex weights)

        # ---------------------------------------------------------
        # 4. Solve the Pricing Problem (Maximum Weight Independent Set)
        # ---------------------------------------------------------
        # Maximize sum(pi_v * y_v) -> Minimize -sum(pi_v * y_v)
        c_pricing = -pi
        res_pricing = milp(c=c_pricing, constraints=pricing_constraints,
                           integrality=integrality, bounds=bounds)

        if not res_pricing.success:
            raise RuntimeError("Pricing MILP failed to converge.")

        W = -res_pricing.fun # Maximum weight of any independent set
        y = np.round(res_pricing.x).astype(int) # The actual independent set

        # ---------------------------------------------------------
        # 5. Compute the Valid Lower Bound
        # ---------------------------------------------------------
        current_lb = Z_rmp / W if W > 0 else Z_rmp
        best_lower_bound = max(best_lower_bound, current_lb)

        print(f"{iteration:<5} | {Z_rmp:<15.4f} | {W:<15.4f} | {best_lower_bound:<15.4f}")

        # ---------------------------------------------------------
        # 6. Check for Optimality
        # ---------------------------------------------------------
        # If no independent set has weight > 1, we have found the exact solution.
        if W <= 1.000001:
            print("\nOptimal fractional chromatic number found!")
            return Z_rmp, True

        # Add the new independent set as a column to the Master matrix
        new_col = y.reshape(-1, 1)
        A_master = np.hstack([A_master, new_col])

    print("\nReached maximum iterations.")
    return best_lower_bound, False

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import time

def greedy_mwis(G, weights):
    """
    Fast greedy heuristic to find an independent set with weight > 1.
    Returns the set as a bit-vector and its total weight.
    """
    nodes = sorted(G.nodes(), key=lambda v: weights[v] / (G.degree(v) + 1), reverse=True)
    ind_set = set()
    total_weight = 0
    remaining_nodes = set(G.nodes())

    for v in nodes:
        if v in remaining_nodes:
            ind_set.add(v)
            total_weight += weights[v]
            # Remove node and all its neighbors
            remaining_nodes.remove(v)
            for neighbor in G.neighbors(v):
                remaining_nodes.discard(neighbor)

    # Convert to vector
    vec = np.zeros(len(G.nodes()))
    for v in ind_set:
        vec[v] = 1
    return vec, total_weight

def fast_fractional_chromatic_bound(G, max_iter=200, time_limit=300):
    n = G.number_of_nodes()
    nodes = list(G.nodes())

    # Pre-build pricing constraints for MILP (the fallback)
    edges = list(G.edges())
    row, col, data = [], [], []
    for i, (u, v) in enumerate(edges):
        row.extend([i, i]); col.extend([u, v]); data.extend([1, 1])
    A_pricing = sp.csr_matrix((data, (row, col)), shape=(len(edges), n))
    pricing_constraints = LinearConstraint(A_pricing, ub=np.ones(len(edges)))

    # Master Problem Init
    A_master = np.eye(n)
    best_lb = 0.0
    start_time = time.time()

    print(f"Accelerated Column Generation: {n} nodes, {len(edges)} edges")
    print(f"{'Iter':<5} | {'Z_RMP':<10} | {'Method':<10} | {'W':<10} | {'Lower Bound':<10}")

    for iteration in range(max_iter):
        if time.time() - start_time > time_limit: break

        # 1. Solve Master Dual
        res_dual = linprog(-np.ones(n), A_ub=A_master.T, b_ub=np.ones(A_master.shape[1]),
                           bounds=(0, None), method='highs')
        if not res_dual.success: break

        Z_rmp = -res_dual.fun
        pi = res_dual.x

        # 2. Try Greedy Heuristic first (FAST)
        y, W = greedy_mwis(G, pi)
        method = "Greedy"

        # 3. Fallback to MILP only if Greedy fails to find a "violating" column
        if W <= 1.001:
            res_milp = milp(c=-pi, constraints=pricing_constraints,
                            integrality=np.ones(n), bounds=Bounds(0, 1))
            if res_milp.success:
                y = np.round(res_milp.x)
                W = -res_milp.fun
                method = "MILP"

        # 4. Calculate Valid Lower Bound
        # ONLY update best_lb if we used the MILP solver (exact pricing)
        if method == "MILP":
            current_lb = Z_rmp / W if W > 1 else Z_rmp
            best_lb = max(best_lb, current_lb)
        else:
            # If greedy, we don't have a mathematically proven lower bound yet
            current_lb = float('nan')

        print(f"{iteration:<5} | {Z_rmp:<10.3f} | {method:<10} | {W:<10.3f} | {best_lb:<10.3f}")

        # 5. Check for convergence
        if W <= 1.0000001:
            print("\nOptimal fractional chromatic number found!")
            best_lb = Z_rmp
            is_exact = True
            break

        A_master = np.hstack([A_master, y.reshape(-1, 1)])

    return best_lb, is_exact

# ==========================================
# Example: Uniform Random Graph G(n, 0.1)
# ==========================================
if __name__ == "__main__":
    # Create G(500, 0.1)
    # Seed fixed for reproducible demonstration
    G = nx.erdos_renyi_graph(100, 0.1, seed=42)

#    lb, is_exact = compute_fractional_chromatic_bound(G, max_iter=2000, time_limit=3600)
    lb, is_exact = fast_fractional_chromatic_bound(G, max_iter=100000, time_limit=3600)

    print("-" * 65)
    if is_exact:
        print(f"EXACT Fractional Chromatic Number: {lb:.4f}")
    else:
        print(f"BEST LOWER BOUND achieved: {lb:.4f}")
