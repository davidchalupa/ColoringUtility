import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import time

# ----------------------
# Preliminary algorithms
# ----------------------

def fractional_chromatic_column_generation_pure(G, max_iter=100000, time_limit=360000):
    """
    Computes a tight lower bound for the fractional chromatic number
    using Column Generation.
    """
    n = G.number_of_nodes()
    if n == 0:
        return 0.0
    start_time = time.time()

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
        elapsed = time.time() - start_time
        if elapsed > time_limit:
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
        if W <= 1.0000001:
            print("-" * 85)
            print(f"\n✅ PROVEN OPTIMALITY REACHED")
            print(f"Exact Fractional Chromatic Number: {Z_rmp:.4f}")
            print(f"Total Time: {elapsed:.2f} seconds")
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

def fractional_chromatic_column_generation_with_greedy(G, max_iter=100000, time_limit=360000):
    n = G.number_of_nodes()
    start_time = time.time()

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
        elapsed = time.time() - start_time

        if time_limit and elapsed > time_limit:
            print("\nTime limit reached.")
            break

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

        print(f"{iteration:<5} | {Z_rmp:<10.4f} | {method:<10} | {W:<10.4f} | {best_lb:<10.4f}")

        # 5. Check for convergence
        if W <= 1.0000001:
            print("-" * 85)
            print(f"\n✅ PROVEN OPTIMALITY REACHED")
            print(f"Exact Fractional Chromatic Number: {Z_rmp:.4f}")
            print(f"Total Time: {elapsed:.2f} seconds")
            return Z_rmp, True

        A_master = np.hstack([A_master, y.reshape(-1, 1)])

    return best_lb, is_exact

# ------------------
# Advanced algorithm
# ------------------

# ==========================================
# 1. HEURISTICS & PRICING
# ==========================================

def dsatur_warm_start(G):
    """Generates an initial feasible set of columns using DSATUR coloring."""
    n = G.number_of_nodes()
    coloring = nx.greedy_color(G, strategy='DSATUR')
    cols = []

    for c_idx in range(max(coloring.values()) + 1):
        v = np.zeros(n)
        for node, color in coloring.items():
            if color == c_idx: v[node] = 1
        cols.append(v)

    # Ensure full mathematical feasibility by adding singletons
    for i in range(n):
        v = np.zeros(n); v[i] = 1; cols.append(v)
    return cols

def greedy_mwis_batch(G, weights, num_trials=20):
    """Tier 1: Super-fast randomized greedy search. Returns multiple columns."""
    n = G.number_of_nodes()
    found_cols = []

    for _ in range(num_trials):
        # Add slight noise to weights to explore different paths
        w_noisy = weights * np.random.uniform(0.9, 1.1, n)
        order = np.argsort(-w_noisy)

        curr_set = np.zeros(n, dtype=bool)
        for v in order:
            if not any(curr_set[nb] for nb in G.neighbors(v)):
                curr_set[v] = True

        # If valid for Master problem
        if np.dot(curr_set, weights) > 1.0001:
            found_cols.append(curr_set.astype(float))

    # Deduplicate
    unique_cols = list({v.tobytes(): v for v in found_cols}.values())
    max_w = max([np.dot(c, weights) for c in unique_cols]) if unique_cols else 0
    return unique_cols, max_w

def tabu_mwis(G, weights, max_steps=150):
    """Tier 2: Tabu search. Used when greedy gets stuck in local optima."""
    n = len(weights)
    order = np.argsort(-weights)

    current_set = np.zeros(n, dtype=bool)
    for v in order:
        if not any(current_set[nb] for nb in G.neighbors(v)):
            current_set[v] = True

    best_set = current_set.copy()
    best_w = np.dot(best_set, weights)
    tabu_list = []

    for _ in range(max_steps):
        improved = False
        candidates = [v for v in range(n) if v not in tabu_list]
        np.random.shuffle(candidates)

        for v in candidates[:40]: # Check a subset to remain fast
            neighbors_in = [nb for nb in G.neighbors(v) if current_set[nb]]
            gain = weights[v] - sum(weights[nb] for nb in neighbors_in)

            if gain > 1e-5:
                current_set[v] = True
                for nb in neighbors_in: current_set[nb] = False

                tabu_list.append(v)
                if len(tabu_list) > 15: tabu_list.pop(0)

                curr_w = np.dot(current_set, weights)
                if curr_w > best_w:
                    best_w = curr_w
                    best_set = current_set.copy()
                improved = True
                break

        if not improved: break

    return best_set.astype(float), best_w

# ==========================================
# 2. MAIN ALGORITHM
# ==========================================

def fractional_chromatic_column_generation_fast(G, max_iter=None, time_limit=None):
    n = G.number_of_nodes()
    start_time = time.time()

    # Pre-build MILP pricing constraints
    edges = list(G.edges())
    row, col, data = [], [], []
    for i, (u, v) in enumerate(edges):
        row.extend([i, i]); col.extend([u, v]); data.extend([1, 1])
    A_pricing = sp.csr_matrix((data, (row, col)), shape=(len(edges), n))
    pricing_constraints = LinearConstraint(A_pricing, ub=np.ones(len(edges)))
    milp_integrality = np.ones(n)
    milp_bounds = Bounds(0, 1)

    # Initialize Master Pool
    cols = dsatur_warm_start(G)

    best_lower_bound = 0.0
    pi_stable = np.zeros(n)
    alpha = 0.5 # Smoothing factor
    iteration = 0

    print(f"Graph: {n} nodes, {len(edges)} edges, Density: {nx.density(G):.3f}")
    print("-" * 85)
    print(f"{'Iter':<5} | {'Z_RMP (Upper)':<14} | {'Method':<7} | {'W_max':<8} | {'Lower Bound':<12} | {'Time'}")
    print("-" * 85)

    while True:
        iteration += 1
        elapsed = time.time() - start_time

        if max_iter and iteration > max_iter:
            print("\nMax iterations reached.")
            break
        if time_limit and elapsed > time_limit:
            print("\nTime limit reached.")
            break

        # 1. Solve Master Problem (Dual Formulation)
        A_dual = np.array(cols)
        res_lp = linprog(-np.ones(n), A_ub=A_dual, b_ub=np.ones(len(cols)),
                         bounds=(0, None), method='highs-ds')

        if not res_lp.success:
            print("\nMaster LP failed to converge.")
            break

        Z_rmp = -res_lp.fun
        pi_raw = res_lp.x

        # 2. Dual Stabilization
        pi_stable = alpha * pi_stable + (1 - alpha) * pi_raw

        # 3. Tiered Pricing Strategy
        new_cols = []

        # Tier 1: Try Fast Greedy
        greedy_cols, W = greedy_mwis_batch(G, pi_stable)
        method = "Greedy"
        if W > 1.0001:
            new_cols.extend(greedy_cols)
        else:
            # Tier 2: Try Tabu Search
            y_tabu, W = tabu_mwis(G, pi_stable)
            method = "Tabu"
            if W > 1.0001:
                new_cols.append(y_tabu)
            else:
                # Tier 3: The Heavy Exact Solver (MILP)
                # Only run when heuristics tap out. Uses raw pi for mathematical proof.
                method = "Exact"
                res_milp = milp(c=-pi_raw, constraints=pricing_constraints,
                                integrality=milp_integrality, bounds=milp_bounds)

                if res_milp.success:
                    y_milp = np.round(res_milp.x)
                    W = -res_milp.fun

                    # Update rigorous mathematical lower bound
                    current_lb = Z_rmp / max(1.0, W)
                    best_lower_bound = max(best_lower_bound, current_lb)

                    if W > 1.0001:
                        new_cols.append(y_milp)
                else:
                    print("MILP solver failed.")
                    break

        # 4. Logging
        # We display the actual lower bound whenever the exact solver proves it.
        # Otherwise, we calculate a temporary "Lagrangian" estimation based on the heuristic W.
        display_lb = best_lower_bound if method == "Exact" else Z_rmp / max(1.0, W)

        if iteration % 5 == 0 or method == "Exact":
            print(f"{iteration:<5} | {Z_rmp:<14.4f} | {method:<7} | {W:<8.4f} | {display_lb:<12.4f} | {elapsed:5.1f}s")

        # 5. Check for Optimality
        if method == "Exact" and W <= 1.00001:
            print("-" * 85)
            print(f"\n✅ PROVEN OPTIMALITY REACHED")
            print(f"Exact Fractional Chromatic Number: {Z_rmp:.4f}")
            print(f"Total Time: {elapsed:.2f} seconds")
            return Z_rmp, True

        # Add new columns to the pool
        cols.extend(new_cols)

    # If it breaks out early via limits
    print("-" * 85)
    print(f"\n⚠️ HALTED EARLY")
    print(f"Best Valid Lower Bound: {best_lower_bound:.4f}")
    print(f"Current Upper Bound (Z_RMP): {Z_rmp:.4f}")
    return best_lower_bound, False

# ==========================================
# 3. EXECUTION
# ==========================================
if __name__ == "__main__":
    G = nx.erdos_renyi_graph(500, 0.1, seed=42)

#    n = 100
#    w = 4
#    G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#    final_val, is_optimal = fractional_chromatic_column_generation_pure(G)
#    final_val, is_optimal = fractional_chromatic_column_generation_with_greedy(G)
    final_val, is_optimal = fractional_chromatic_column_generation_fast(G, max_iter=None, time_limit=None)
