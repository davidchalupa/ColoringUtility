import cvxpy as cp
import networkx as nx
import numpy as np
import time

def lovasz_theta(G):
    n = G.number_of_nodes()

    print(f"[*] Phase 1: Preparing Adjacency and Mask Matrices ({n} vertices)...")
    start_time = time.time()

    # We calculate the number of complement edges mathematically for the diagnostic output
    num_edges_comp = int((n * (n - 1) / 2) - G.number_of_edges())
    print(f"    - Implicit complement edges (masked): {num_edges_comp}")

    # Create a Mask for the complement edges
    # adj is the adjacency matrix of G.
    # mask is 1 where edges exist in G, plus the diagonal (identity matrix).
    adj = nx.to_numpy_array(G)
    mask = adj + np.eye(n)

    # Define the SDP variable
    X = cp.Variable((n, n), PSD=True)

    print("[*] Phase 2: Constructing Vectorized SDP Constraints...")
    # The vectorized constraint replaces the slow for-loop
    # cp.multiply does element-wise multiplication. (1 - mask) is 1 only where
    # edges are MISSING in G, forcing those specific X[i,j] values to 0.
    constraints = [
        cp.trace(X) == 1,
        cp.multiply(X, 1 - mask) == 0
    ]

    # Objective: Maximize the sum of all entries in X
    prob = cp.Problem(cp.Maximize(cp.sum(X)), constraints)

    print("[*] Phase 3: Handing over to Solver (SCS)...")
    print("-" * 30)
    # verbose=True will print the iteration-by-iteration progress of SCS
    # eps=1e-6 for higher precision
    prob.solve(solver=cp.SCS, verbose=True, eps=1e-6)
    print("-" * 30)

    end_time = time.time()
    print(f"[*] Done! Total time: {end_time - start_time:.2f} seconds")

    return prob.value

G = nx.erdos_renyi_graph(250, 0.1, seed=42)
lb = lovasz_theta(G)
print(f"Lovász Lower Bound: {lb:.4f}")
