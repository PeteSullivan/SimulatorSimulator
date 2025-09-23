import GateGeneration as gg
import MathematicaScripts as scripts
import numpy as np
from scipy.optimize import minimize, basinhopping
from scipy.optimize import differential_evolution, minimize

def optimize_frequencies_direct(Jij_ideal, equilibrium_positions, secular_frequencies, NIons, n_frequencies, 
                               freq_bounds=(3.0, 7.0), max_iter=1000, tol=1e-15):
    """
    Directly optimizes multiple laser frequencies to minimize infidelity with the ideal interaction.
    
    Args:
        equilibrium_positions: Ion equilibrium positions
        secular_frequencies: Secular frequencies of the trap
        NIons: Number of ions
        n_frequencies: Number of frequencies to optimize
        freq_bounds: Frequency search range in MHz
        max_iter: Maximum optimization iterations
        tol: Convergence tolerance
    
    Returns:
        tuple: (optimal_frequencies, optimal_amplitudes, best_infidelity)
    """
    # Get mode frequencies in MHz for initialization
    mode_freqs, mode_vectors = gg.eigensystem(equilibrium_positions, secular_frequencies, NIons)
    
    # Initialize frequencies spread across the range
    freq_min, freq_max = freq_bounds
    initial_freqs = np.linspace(freq_min, freq_max, n_frequencies)
    
    # Use same initial amplitude for all lasers
    initial_amplitude = 180e3 / np.sqrt(n_frequencies)  
    initial_amplitudes = np.linspace(initial_amplitude, 180e3, n_frequencies)
    
    # Combined parameters vector
    initial_params = np.concatenate([initial_freqs, initial_amplitudes])
    
    # Define bounds
    bounds = [(freq_min, freq_max)] * n_frequencies + [(10e3, 500e3)] * n_frequencies
    
    def objective(params):
        """Objective function that directly calculates infidelity"""
        freqs = params[:n_frequencies]
        amplitudes = params[n_frequencies:]
        
        # Create lasers
        lasers = []
        for i, (freq, amp) in enumerate(zip(freqs, amplitudes)):
            laser = gg.Laser(
                frequency=freq * 1e6,  # Convert to Hz
                Omegas=[amp] * NIons,
                phase=0
            )
            lasers.append(laser)

        # Get the resulting interaction matrix
        Jij_exp = gg.Jij(
            mode_freqs,
            mode_vectors,
            NIons,
            lasers
        )

        # Calculate infidelity with the ideal interaction
        inf = Infidelity(Jij_exp, Jij_ideal)
        
        return inf
    
    # Run optimization
    result = minimize(objective, initial_params, bounds=bounds, method='TNC', 
                     options={'maxiter': max_iter, 'ftol': tol})

    # result = basinhopping(objective, initial_params, target_accept_rate=0.001)
    
    optimal_freqs = result.x[:n_frequencies]
    optimal_amplitudes = result.x[n_frequencies:]
    best_infidelity = result.fun
    
    return optimal_freqs, optimal_amplitudes, best_infidelity

def Jij_Generation(equilibrium_positions, secular_frequencies, NIons, Jij_ideal):
    """
    Optimizes 2*NIons+1 laser frequencies to generate the desired interaction pattern.
    Uses direct infidelity minimization as the objective function.
    """

    # Number of frequencies to optimize
    n_frequencies = 2 * NIons + 1
    
    # Get mode frequencies in MHz for setting bounds
    mode_freqs, _ = gg.eigensystem(equilibrium_positions, secular_frequencies, NIons)
    mode_freqs_MHz = mode_freqs / 1e6
    print("frequencies:", mode_freqs_MHz)
    # Find optimal frequencies and amplitudes
    opt_freqs, opt_amplitudes, infidelity = optimize_frequencies_direct(
        Jij_ideal,
        equilibrium_positions,
        secular_frequencies,
        NIons,
        n_frequencies=n_frequencies,
        freq_bounds=(min(mode_freqs_MHz)-0.1, max(mode_freqs_MHz)+0.1)
    )
    # Create lasers with the optimized parameters
    lasers = []
    print(f"Optimized {n_frequencies} frequencies:")
    for i, (freq, amp) in enumerate(zip(opt_freqs, opt_amplitudes)):
        laser = gg.Laser(
            frequency=freq * 1e6,  # Convert to Hz
            Omegas=[amp] * NIons,
            phase=0
        )
        lasers.append(laser)
        print(f"  Laser {i+1}: {freq:.3f} MHz, amplitude: {amp/1000:.1f} kHz")
    
    print(f"Final infidelity: {infidelity:.6f}")
    
    return lasers, infidelity


def GlobalInteractionIdeal(NIons):
    Jij = np.zeros([NIons, NIons])
    for i in range(NIons):
        for j in range(NIons):
            if i != j:
                Jij[i][j] = 1
    return Jij

def Infidelity(J_exp, J_des):
    """
    Calculates the infidelity between two J matrices using the formula:

        I ≡ 1/2 * (1 - ⟨J_exp, J_des⟩ / (∥J_exp∥ ∥J_des∥))

    Parameters:
    J_exp (array): Experimental J matrix.
    J_des (array): Desired J matrix.

    Returns:
    float: Infidelity value.
    """
    # Remove diagonal elements
    J_exp_off_diag = J_exp - np.diag(np.diag(J_exp))
    J_des_off_diag = J_des - np.diag(np.diag(J_des))

    # Compute the Frobenius inner product of the two matrices
    frobenius_product = np.sum(J_exp_off_diag * J_des_off_diag)

    # Compute the Frobenius norms of the matrices
    norm_exp = np.sqrt(np.sum(J_exp_off_diag * J_exp_off_diag))
    norm_des = np.sqrt(np.sum(J_des_off_diag * J_des_off_diag))

    # Compute the infidelity
    infidelity_value = 0.5 * (1 - (frobenius_product / (norm_exp * norm_des)))
    return infidelity_value

# def Infidelity(J_exp, J_des, alpha=0.01, epsilon=1e-10):
#     """
#     Enhanced infidelity metric with normalization, weighting, and regularization.
    
#     Args:
#         J_exp (np.ndarray): Experimental J matrix
#         J_des (np.ndarray): Desired J matrix
#         alpha (float): Regularization strength (default: 0.01)
#         epsilon (float): Small constant to avoid division by zero
        
#     Returns:
#         float: Composite infidelity measure (0 = perfect match)
#     """
#     # Remove diagonal elements (self-interactions typically not important)
#     J_exp_off = J_exp - np.diag(np.diag(J_exp))
#     J_des_off = J_des - np.diag(np.diag(J_des))
    
#     # Normalization factors
#     norm_exp = np.linalg.norm(J_exp_off, 'fro')
#     norm_des = np.linalg.norm(J_des_off, 'fro')
    
#     # Core infidelity metric (normalized inner product)
#     core_infidelity = 0.5 * (1 - np.sum(J_exp_off * J_des_off) / 
#                      (norm_exp * norm_des + epsilon))
    
#     # Relative error term (emphasizes large relative errors)
#     relative_error = np.mean(np.abs(J_exp_off - J_des_off) / 
#                            (np.abs(J_des_off) + epsilon))
    
#     # Structure preservation term (maintains sign patterns)
#     sign_agreement = np.mean(np.sign(J_exp_off) != np.sign(J_des_off))
    
#     # Regularization term (penalizes extreme amplitudes)
#     regularization = alpha * (np.mean(np.abs(J_exp_off)) / norm_des)
    
#     # Composite infidelity measure
#     composite_infidelity = (
#         core_infidelity + 
#         0.3 * relative_error + 
#         0.2 * sign_agreement + 
#         regularization
#     )
    
#     return composite_infidelity



#Jij_Generation_NEW pseudocode
#1. get frequencies for beams (assume correct for now, needs fixing?)
#2. build J_k matrices



def Jij_Generation_NEW(equilibrium_positions, secular_frequencies, NIons, Jij_ideal):
    """
    Optimizes 2*NIons+1 laser frequencies to generate the desired interaction pattern.
    Uses direct infidelity minimization as the objective function.
    """
    #assuming we have the correct frequencies, only optimize the amplitudes
    #each frequency affects each mode with a certain amount, and can only be scaled proportionally
    #we need to find how much of each laser frequency needed to get the right proportions of each mode
    #Ax = b, A is the frequencies 

    
    # Number of frequencies to optimize
    n_frequencies = 2 * NIons + 1
    
    
    
    # Get mode frequencies in MHz for setting bounds
    mode_freqs, mode_vectors = gg.eigensystem(equilibrium_positions, secular_frequencies, NIons)
    mode_freqs_MHz = mode_freqs / 1e6
    print("frequencies:", mode_freqs_MHz)
    # Find optimal frequencies and amplitudes
    opt_freqs, opt_amplitudes, infidelity = optimize_frequencies_direct(
        Jij_ideal,
        equilibrium_positions,
        secular_frequencies,
        NIons,
        n_frequencies=n_frequencies,
        freq_bounds=(min(mode_freqs_MHz)-0.1, max(mode_freqs_MHz)+0.1)
    )
    print("opt_freqs:", opt_freqs)
    
    #for multimode, less important to do rn
    # #build J_k matrices from flipped qubits
    # J_ks = []
    # for positive in range(2): #second half of weights are negative
    #     for k in range(NIons-1): #different mode vectors
    #         for n in range(NIons+1): #flipped qubit
    #             J_k = np.zeros([NIons, NIons])
    #             for i in range(NIons):
    #                 for j in range(NIons):
    #                     J_k[i][j] = mode_vectors[k][i] * mode_vectors[k][j]
    #                     if i == n-1 or j == n-1:
    #                         J_k[i][j] = -J_k[i][j]
    #             if positive == 1:
    #                 J_k = -J_k
    #             # print(f"positive? {positive}, mode vector: {k}, flipped: {n}")
    #             # gg.print_matrix(J_k)
    #             J_ks.append(J_k)

    #for single mode
    J_ks = []
    for i in range(NIons):
        J_k = np.zeros((NIons, NIons))
        for m in range(NIons):
            for n in range(NIons):
                J_k[m][n] = mode_vectors[i][m] * mode_vectors[i][n]
        J_ks.append(J_k)
    


    # FIND WEIGHTS FOR DIFFERENT J_k MATRIX:    
    # Flatten each matrix into a vector
    A = np.column_stack([M.flatten() for M in J_ks])  # shape: (m*k, n)
    b = Jij_ideal.flatten()  # shape: (m*k,)
    
    # Solve linear system A w = b
    weights, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    
    print("weights:", weights)

    #make sure weights give correct Jij
    #multiply each by the weights
    Jij = np.zeros([NIons, NIons])
    for i in range(NIons):
        Jij += J_ks[i] * weights[i]
    #output Jij
    print("Jij generated:", Jij)
    print("Jij ideal:", Jij_ideal)
    # return Jij




    #create interaction matrix between opt_freqs and modes, A_ij = ith mode affected by jth frequency
    Interaction_Strengths = np.zeros((NIons, len(opt_freqs)))
    for i in range(NIons):
        for j in range(len(opt_freqs)):
            Interaction_Strengths[i][j] = 1/(opt_freqs[j]**2-mode_freqs_MHz[i])
    print("Interaction_Strengths matrix:")
    gg.print_matrix(Interaction_Strengths)

    x, residuals, rank, s = np.linalg.lstsq(Interaction_Strengths, weights, rcond=None)
    norm = np.linalg.norm(x)
    normalized_weights = x / norm
    print("x:", normalized_weights)
    
    # Create lasers with the optimized parameters
    lasers = []
    print(f"Optimized {n_frequencies} frequencies:")
    for i, (freq, amp) in enumerate(zip(opt_freqs, x)):
        laser = gg.Laser(
            frequency=freq * 1e6,  # Convert to Hz
            Omegas=[amp] * NIons,
            phase=0
        )
        lasers.append(laser)
        print(f"  Laser {i+1}: {freq:.3f} MHz, amplitude: {amp/1000:.1f} kHz")
    
    print(f"Final infidelity: {infidelity:.6f}")
    
    return lasers, infidelity