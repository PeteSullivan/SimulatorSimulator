import numpy as np
from scipy.optimize import minimize

def find_optimal_weights(J_desired, J_modes):
    """
    Find the optimal weights c_k to approximate J_desired using the mode matrices J_modes.
    
    Parameters:
        J_desired (np.ndarray): Desired spin-spin coupling matrix (N x N).
        J_modes (list[np.ndarray]): List of mode interaction matrices J^(k) (each N x N).
    
    Returns:
        np.ndarray: Optimal weights c_k.
        float: Infidelity (1 - fidelity) between J_desired and the approximation.
    """
    N = J_desired.shape[0]
    num_modes = len(J_modes)
    
    # Flatten J_desired and J_modes for optimization
    J_des_flat = J_desired.flatten()
    J_modes_flat = np.array([Jk.flatten() for Jk in J_modes])
    
    # Objective function: minimize infidelity (Frobenius norm difference)
    def objective(c):
        J_approx = np.sum(c[:, None] * J_modes_flat, axis=0)
        diff = J_approx - J_des_flat
        return np.sum(diff**2)  # Frobenius norm squared
    
    # Constraint: sum(c_k) = 0 (to avoid identity matrix contribution)
    constraints = {'type': 'eq', 'fun': lambda c: np.sum(c)}
    
    # Bounds: c_k can be positive or negative
    bounds = [(-np.inf, np.inf) for _ in range(num_modes)]
    
    # Initial guess: random weights
    c0 = np.random.randn(num_modes)
    
    # Optimize
    result = minimize(objective, c0, bounds=bounds, constraints=constraints)
    c_opt = result.x
    
    # Compute infidelity (Eq. 12 in the paper)
    J_approx = np.sum(c_opt[:, None, None] * J_modes, axis=0)
    norm_J_des = np.linalg.norm(J_desired, 'fro')
    norm_J_approx = np.linalg.norm(J_approx, 'fro')
    inner_prod = np.sum(J_desired * J_approx)
    
    infidelity = 0.5 * (1 - inner_prod / (norm_J_des * norm_J_approx))
    
    return c_opt, infidelity

# Example usage:
if __name__ == "__main__":
    # Example: 3 ions, 2 modes (COM and stretch)
    N = 3
    J_modes = [
        np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]),  # COM mode (all-to-all)
        np.array([[1, 0, -1], [0, 0, 0], [-1, 0, 1]])  # Stretch mode
    ]
    
    # Desired J matrix (e.g., nearest-neighbor couplings)
    J_desired = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])
    
    c_opt, infidelity = find_optimal_weights(J_desired, J_modes)
    print("Optimal weights c_k:", c_opt)
    print("Infidelity:", infidelity)