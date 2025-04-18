import numpy as np
import math
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint
# from scipy.optimize import minimize, curve_fit, fsolve
# import scipy.constants as scipy
import time
# import qutip as qt

class Laser:
    '''
    a laser in the system used to create gates. omegas correspond to each ion's rabi frequency.
    '''
    def __init__(self, frequency, Omegas=None, phase=0, n=None, c=0, max_Omega=1.0, sigma=1.0):
        self.frequency = frequency
        self.phase = phase

        if Omegas is not None:
            # Use predefined Omegas if provided
            self.Omegas = Omegas
        else:
            # Generate Gaussian-distributed Omegas if no list is given
            if n is None:
                raise ValueError("If Omegas is not provided, 'n' (number of ions) must be specified.")
            self.Omegas = self._generate_gaussian_omegas(n, c, max_Omega, sigma)

    def _generate_gaussian_omegas(self, n, c, max_Omega, sigma):
        '''
        Generates a Gaussian distribution of Rabi frequencies.
        
        Args:
            n (int): Number of ions.
            c (int): Center index (peak of Gaussian).
            max_Omega (float): Maximum Rabi frequency at center.
            sigma (float): Standard deviation of the Gaussian.
        
        Returns:
            list: Gaussian-distributed Rabi frequencies.
        '''
        omegas = []
        for i in range(n):
            # Gaussian function: max_Omega * exp(-(i - c)^2 / (2 sigma^2))
            omega = max_Omega * math.exp(-((i - c) ** 2) / (2 * sigma ** 2))
            omegas.append(omega)
        return omegas



# region MS gate

def find_directional_frequencies(equilibrium_positions, secular_frequencies):

    '''
    finds the transverse direction for 2d lattices, used to find the eigensystem
    of a trap potential.
    '''
    #finds which direction is all 0, return the correct frad and fax to use
    #also converts MHz to Hz for eigensystem function
    for i in range(len(equilibrium_positions)):
        if (equilibrium_positions[i][0] != 0):
            return secular_frequencies[1]*10**6, secular_frequencies[2]*10**6
    return secular_frequencies[0], secular_frequencies[2]

def eigensystem(equilibrium_positions, secular_frequencies, NIons):
    
    """
    Calculates the mode spectrum for a given set of ions at specific trap frequencies.

    Parameters:
    frad (float): Radial Trap frequency in MHz.
    fax (float): Axial trap frequency in MHz.
    num_ions (int): The number of ions in the system.
    lc (float): The lattice constant in meters.

    Returns:
    eigenvalues (array): Normal mode frequencies in MHz. Need to be
    multiplied by 10^6 to get in Hz.

    eigenvectors (array): Array of eigenvectors.
    """
    e = 1.602176634e-19  # Elementary charge in C
    eps_o = 8.854187817e-12  # Vacuum permittivity in F/m
    m = 171 * 1.67262192e-27  # Mass of Yb-171 ion in kg
    frad, fax = find_directional_frequencies(equilibrium_positions, secular_frequencies)
    lc = (e**2 / (4 * np.pi * eps_o * m * fax**2)) ** (1/3)
    u = equilibrium_positions

    Anm = np.empty((NIons,NIons))
    # print("u:", u)
    for n in range(NIons):
        for m in range(NIons):
            if n == m:
                sum1 = 0.0
                for p in range(NIons):
                    if(p!=m):
                        # sum1 += 1 / abs(np.sqrt((u[m][0]/lc - u[p][0]/lc)**2 + (u[m][1]/lc - u[p][1]/lc)**2 + (u[m][2]/lc - u[p][2]/lc)**2))**3
                        sum1 += 1 / abs(np.sqrt((u[m][0] - u[p][0])**2 + (u[m][1] - u[p][1])**2 + (u[m][2] - u[p][2])**2)/lc)**3
                Anm[n][n]= (frad/fax)**2-sum1
                # print("frad/fax:", frad/fax)
            elif n!=m:
                # sum2 = 1 / abs(np.sqrt((u[m][0]/lc - u[n][0]/lc)**2 + (u[m][1]/lc - u[n][1]/lc)**2 + (u[m][2]/lc - u[n][2]/lc)**2))**3
                sum2 = 1 / abs(np.sqrt((u[m][0] - u[n][0])**2 + (u[m][1] - u[n][1])**2 + (u[m][2] - u[n][2])**2)/lc)**3
                Anm[n][m]= sum2
    # print("Anm:", Anm)
    eigenvalues, eigenvectors = np.linalg.eig(Anm)
    eigenvalues = np.sqrt(eigenvalues)*fax
    eigenvectors = eigenvectors.T
    eigenvectors = eigenvectors / np.linalg.norm(eigenvectors, axis=0, keepdims=True)

    #sort the modes by eigenvalues
    sorted_indices = np.argsort(eigenvalues)[::-1] # COM mode = first mode
    mode_frequencies = np.sort(eigenvalues)[::-1]
    mode_vectors = eigenvectors[sorted_indices]
    # print("mode_frequencies:", mode_frequencies)
    # print("eigenvectors:", eigenvectors)
    return mode_frequencies, mode_vectors

def find_Jks(mode_frequencies, mode_vectors):
    '''
    returns the J(k) matrices for global beam multimode laser interaction engineering
    '''
    return [np.outer(mode_vectors[i, :], mode_vectors[i, ]) for i in range(len(mode_frequencies))]
    
def GlobalBeamJij(mode_frequencies, mode_vectors, NIons, laserset):
    '''
    Calculates the effective Jij matrix from a set of lasers and a trap's mode data
    '''
    
    hbar = 1.054571817e-34  # Reduced Planck's constant (J·s)
    delta_k = np.sqrt(2) * 2 * np.pi / (355e-9)
    m = 171 * 1.67262192e-27  # Mass of Yb-171 ion in kg
    recoil = hbar * delta_k**2 / (2 * m)

    # Initialize the total J matrix
    total_J_int = np.zeros((NIons, NIons))

    #find Jk matrices
    Jks = find_Jks(mode_frequencies, mode_vectors)

    #sum all laser effects for each mode to make Jij matrix
    for mode in range(len(mode_frequencies)):
        mode_freq = mode_frequencies[mode]
        for laser in range(len(laserset)):

            c_k =  laserset[laser].Omegas[0] * laserset[laser].Omegas[0] * recoil / ((laserset[laser].frequency)**2 - mode_freq**2)
            total_J_int += Jks[mode] * c_k/ (2*np.pi*10**3)  # Sum all interaction matrices (output in Hz)

    #remove diagonal
    for n in range(NIons):
        total_J_int[n][n] = 0

    # print("Final Jij Matrix:\n", total_J_int)
    return total_J_int
    return np.round(total_J_int, 3)

def Jij(mode_frequencies, mode_vectors, NIons, laserset):
    """
    Calculates the J matrix that dictate ion-ion coupling. Requires the inputs of
    frequencies to be in Hz.
  
    Note: we could modify for in the future to include position calculations.
  
    Parameters:
  
    Omegas (array): Array of rabi frequencies in Hz.
    mode_vectors (array): Array of eigenvectors.
    mode_frequencies (array): Array of eigenvalues.
    mode_detune (float): The mode being detuned from in Hz. e.g. COM.
    det (float): Detuning from mode in Hz.
    recoil (float): Recoil frequency in Hz.
    num_ions (int): Number of ions.
  
    Returns:
    J (array): The J matrix.
    """
    hbar = 1.054571817e-34  # Reduced Planck's constant (J·s)
    delta_k = np.sqrt(2) * 2 * np.pi / (355e-9)
    m = 171 * 1.67262192e-27  # Mass of Yb-171 ion in kg
    recoil = hbar * delta_k**2 / (2 * m)

    J = np.zeros((NIons,NIons),dtype=float)
    for laser in range(len(laserset)):
        for i in range(NIons):
            for j in range(NIons):
                if i < j:
                    s = sum( (mode_vectors[k][i]*mode_vectors[k][j]) / ((laserset[laser].frequency)**2-mode_frequencies[k]**2) \
                        for k in range(NIons))
                    J[i][j] = recoil * laserset[laser].Omegas[i] * laserset[laser].Omegas[j] * s / (2*np.pi*10**3)
                    # J[i][j] = np.round(J[i][j], decimals=3)

    return J

def Ising_interaction(i,j,N, phase):
    """
    Constructs a two-body XX interaction for an N-particle system.
  
    Parameters:
    i (int): The index of the first particle.
    j (int): The index of the second particle.
    N (int): The total number of particles in the system.
    phase (float): Phase φ used in the transformation.
  
    Returns:
    qt.Qobj: X on i and j, identity everywhere else.
    """
    operators = []
    
    # Identity matrix
    I = np.array([[1, 0], [0, 1]])
    
    # Pauli-X and Pauli-Y matrices
    X = np.array([[0, 1], [1, 0]])  # σx
    Y = np.array([[0, -1j], [1j, 0]])  # σy

    # Compute e^σφ = cos(φ)σx + sin(φ)σy
    sigma_phi = np.cos(phase) * X + np.sin(phase) * Y

    for k in range(N):
        if k == i or k == j:
            operators.append(sigma_phi)
        else:
            operators.append(I)

    # Compute the tensor product
    H = np.array([1])  # Start with scalar 1 (acts as identity in Kronecker product)
    for operator in operators:
        H = np.kron(H, operator)

    return H


def MS_Hamiltonian(J,N, phase):

    """
    Calculates the overall Hamiltonian for an N-particle system.
  
    Parameters:
    J (array): The J matrix.
    N (int): The number of particles in the system.
  
    Returns:
    H (array): The Hamiltonian.
    """
    H = []
    for i in range(N):
      for j in range(N):
        if i < j:
          val = 2*np.pi*(J[i][j] * Ising_interaction(i,j,N, phase))
          H.append(val)
    return sum(H)

def MS_gate_Unitary(equilibrium_positions, secular_frequencies, NIons, laserset, time):
    """Calculates the Unitary gate operation given trap parameters and a set of lasers"""
    mode_frequencies, mode_vectors = eigensystem(equilibrium_positions=equilibrium_positions, secular_frequencies=secular_frequencies, NIons=NIons)
    
    Jij_interactions = Jij(mode_frequencies=mode_frequencies, mode_vectors=mode_vectors, NIons=NIons, laserset=laserset)
    # print("MS Jij:\n", Jij_interactions)
    Hamiltonian = MS_Hamiltonian(J=Jij_interactions, N=NIons, phase=laserset[0].phase)
    # print("MS hamiltonian:\n", Hamiltonian)
    Unitary = unitary_evolution(H=Hamiltonian, t=time)
    # print("Unitary:\n", Unitary)
    return Unitary

def Trap_to_Jij(equilibrium_positions, secular_frequencies, NIons, laserset):
    mode_frequencies, mode_vectors = eigensystem(equilibrium_positions=equilibrium_positions, secular_frequencies=secular_frequencies, NIons=NIons)
    Jij_interactions = Jij(mode_frequencies=mode_frequencies, mode_vectors=mode_vectors, NIons=NIons, laserset=laserset)
    return Jij_interactions

def MS_gate_Hamil(equilibrium_positions, secular_frequencies, NIons, laserset, return_Jij=False):
    """Calculates the Unitary gate operation given trap parameters and a set of lasers"""
    mode_frequencies, mode_vectors = eigensystem(equilibrium_positions=equilibrium_positions, secular_frequencies=secular_frequencies, NIons=NIons)

    Jij_interactions = Jij(mode_frequencies=mode_frequencies, mode_vectors=mode_vectors, NIons=NIons, laserset=laserset)
    # print("MS Jij:\n", Jij_interactions)
    Hamiltonian = MS_Hamiltonian(J=Jij_interactions, N=NIons, phase=laserset[0].phase)
    # print("MS hamiltonian:", Hamiltonian)
    if return_Jij:
        return Hamiltonian, Jij_interactions
    return Hamiltonian

# endregion 

def system_frequency():
    #currently finds system frequency for |0> always, should be a function of the state later
    return 12.6e9

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

def global_z(NIons):
    
    # Pauli Z matrix and identity matrix
    sigma_z = np.array([[1, 0], [0, -1]])
    I = np.array([[1, 0], [0, 1]])
    # Initialize the Hamiltonian as a zero matrix
    H = np.zeros((2**NIons, 2**NIons), dtype=complex)

    # Iterate over each qubit
    for i in range(NIons):
        # Construct the operator for the i-th qubit: I ⊗ I ⊗ ... ⊗ σ_z ⊗ ... ⊗ I
        operator = [I] * NIons  # Start with identity matrices for all qubits
        operator[i] = sigma_z  # Replace the i-th identity with σ_z
        # Compute the tensor product
        sigma_z_i = operator[0]
        for j in range(1, NIons):
            sigma_z_i = np.kron(sigma_z_i, operator[j])
        # Add the contribution to the Hamiltonian
        H += (1 / 2) * sigma_z_i

    return H

def stark_shift(Omegas):
    
    # Pauli Z matrix
    sigma_z = np.array([[1, 0], [0, -1]])
    eye = np.eye(2)
    H = np.zeros([2**len(Omegas), 2**len(Omegas)])

    for i in range(len(Omegas)):
        element = Omegas[i]/2 * sigma_z
        local_rotation = np.array([1])
        for j in range(len(Omegas)):
            if i != j:
                local_rotation = np.kron(local_rotation, eye)
            else:
                local_rotation = np.kron(local_rotation, element)
        H += local_rotation
        


    return H

def unitary_evolution(H, t):
    """
    Computes the unitary evolution matrix U(t) = e^(-i H t)
    
    Parameters:
    H (numpy.ndarray): Hamiltonian matrix (can be complex).
    t (float): Evolution time.

    Returns:
    numpy.ndarray: The unitary evolution matrix U(t).
    """
    hbar = 1.054571817e-34
    # Step 1: Eigen decomposition H = VDV^(-1)
    eigenvalues, eigenvectors = np.linalg.eig(H)  # Use eig() for general matrices
    
    # Step 2: Compute e^(-i * eigenvalues * t) for each eigenvalue
    exp_diag = np.diag(np.exp(-1j * eigenvalues * t))
    
    # Step 3: Reconstruct U = V e^(-iDt) V^(-1)
    U = eigenvectors @ exp_diag @ np.linalg.inv(eigenvectors)
    
    return U

def print_matrix(matrix, precision=3, sci=True):

    """
    Nicely print a complex matrix with formatted real and imaginary parts.
    
    Parameters:
        matrix (np.ndarray): The complex matrix to print.
        precision (int): Number of decimal places to show.
        sci (bool): Whether to use scientific notation.
    """
    rows, cols = matrix.shape
    for i in range(rows):
        row_str = ""
        for j in range(cols):
            z = matrix[i, j]
            if sci:
                fmt = f"{{0.real:.{precision}e}}{{0.imag:+.{precision}e}}j"
            else:
                fmt = f"{{0.real:.{precision}f}}{{0.imag:+.{precision}f}}j"
            row_str += fmt.format(z).rjust(24)
        print(row_str)

