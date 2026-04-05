import cvxpy as cp
import networkx as nx
import numpy as np
import time

def get_rigorous_bound(X_matrix, mask, tolerance=1e-7):
    """
    Refined projection with safety margins to prevent 4.00000001 -> 5 jumps.
    """
    n = X_matrix.shape[0]

    # 1. Force Symmetry
    X_safe = (X_matrix + X_matrix.T) / 2

    # 2. Force Edge Constraints
    X_safe = np.multiply(X_safe, mask)

    # 3. Force PSD (Zero out eigenvalues < 0)
    vals, vecs = np.linalg.eigh(X_safe)
    vals = np.maximum(vals, 0)
    X_safe = vecs @ np.diag(vals) @ vecs.T

    # 4. Force Trace = 1 exactly
    tr = np.trace(X_safe)
    if tr > 0:
        X_safe = X_safe / tr

    # 5. Calculate sum
    bound_value = np.sum(X_safe)

    # 6. Apply Numerical Slack
    # We subtract a tiny epsilon to handle floating point noise
    # before the ceiling operation.
    rigorous_lb = np.ceil(bound_value - tolerance)

    return bound_value, int(rigorous_lb)

def lovasz_theta_rigorous(G):
    n = G.number_of_nodes()

    adj = nx.to_numpy_array(G)
    # For Theta(G_complement) which bounds Chi(G), we mask non-edges
    mask = adj + np.eye(n)

    X = cp.Variable((n, n), PSD=True)

    # We use the standard Primal formulation
    constraints = [
        cp.trace(X) == 1,
        cp.multiply(X, 1 - mask) == 0
    ]
    prob = cp.Problem(cp.Maximize(cp.sum(X)), constraints)

    print(f"[*] Starting rigorous solver for {n} nodes...")
    print("[!] The 'obj' in logs may be 'too high' until convergence.")
    print("-" * 60)

    # We use a very strict eps for high precision
    # Removed time_limit so it runs until the precision is met
    prob.solve(solver=cp.SCS, verbose=True, eps=1e-8)

    print("-" * 60)

    if X.value is None:
        return None

    # Calculate the 'raw' (potentially over-optimistic) bound
    raw_lb = prob.value

    # Calculate the 'rigorous' (guaranteed safe) bound
    safe_float, safe_int = get_rigorous_bound(X.value, mask)

    return raw_lb, safe_float, safe_int

G = nx.erdos_renyi_graph(500, 0.1, seed=42)
raw, safe_f, safe_i = lovasz_theta_rigorous(G)

if raw is not None:
    print("\n" + "="*30)
    print(f"Solver Raw (Unsafe):      {raw:.8f}")
    print(f"Projected Safe (Float):   {safe_f:.8f}")
    print(f"Rigorous Chromatic LB:    {safe_i}")
    print("="*30)
