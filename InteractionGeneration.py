import GateGeneration as gg
import MathematicaScripts as scripts
import numpy as np
from scipy.optimize import minimize, basinhopping

def optimize_frequencies_direct(Jij_ideal, equilibrium_positions, secular_frequencies, NIons, n_frequencies, 
                               freq_bounds=(4.0, 6.0), max_iter=1000, tol=1e-15):
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
    print("frequencies:", mode_freqs)
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

