import numpy as np

def system_frequency():
    #currently finds system frequency for |0> always, should be a function of the state later
    return 12.6e10

def sqr(t, laser):
    """
    Computes the quantum evolution term: 
    Ω/2 * (σ+ * exp(-i * (Δ * t + φ)) + σ- * exp(i * (Δ * t + φ)))

    Parameters:
    Omega (float): Rabi frequency (Ω)
    Delta (float): Detuning (Δ)
    t (float): Time (t)
    phi (float): Phase (φ)
    sigma_plus (complex or matrix): Raising operator (σ+)
    sigma_minus (complex or matrix): Lowering operator (σ-)

    Returns:
    complex or matrix: Result of the equation
    """
    hbar = 1.05457182e-34 #currently not multiplying by hbar to avoid rounding errors, might need to change later
    sigma_plus = np.array([[0, 1], [0, 0]])  # Raising operator
    sigma_minus = np.array([[0, 0], [1, 0]])  # Lowering operator

    detuning = system_frequency() - laser.frequency
    exp_factor = np.exp(-1j * (detuning * t + laser.phase))  # e^(-i(Δt + φ))
    exp_factor_conj = np.exp(1j * (detuning * t + laser.phase))  # e^(i(Δt + φ))
    
    NIons = len(laser.Omegas)
    H = np.zeros((NIons**2,NIons**2),dtype=float)
    print("NIons:", NIons)
    for i in range(NIons):
        H_temp = np.array([1])
        Hamil = (laser.Omegas[i] / 2) * (sigma_plus * exp_factor + sigma_minus * exp_factor_conj)
        for j in range(NIons):
            if i == j:
                H_temp = np.kron(H_temp, Hamil)
            else:
                H_temp = np.kron(H_temp, np.eye(2))
        # print("H:", H)
        # print("H_temp:", H_temp)
        H = H + H_temp
    
    # return H
    return np.round(H, 9)


def unitary_evolution(H, t):
    """
    Computes the unitary evolution matrix U(t) = e^(-i H t)
    
    Parameters:
    H (numpy.ndarray): Hamiltonian matrix (Hermitian).
    t (float): Evolution time.

    Returns:
    numpy.ndarray: The unitary evolution matrix U(t).
    """
    # Step 1: Eigen decomposition H = VDV^(-1)
    eigenvalues, eigenvectors = np.linalg.eigh(H)  # Since H is Hermitian, use eigh()
    
    # Step 2: Compute e^(-i * eigenvalues * t) for each eigenvalue
    exp_diag = np.diag(np.exp(-1j * eigenvalues * t))
    
    # Step 3: Reconstruct U = V e^(-iDt) V^(-1)
    U = eigenvectors @ exp_diag @ np.linalg.inv(eigenvectors)
    
    # return U
    return np.round(U, 9)

