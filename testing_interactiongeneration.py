import MathematicaScripts as scripts
import InteractionGeneration as ig
import GateGeneration as gg
import numpy as np
from scipy.spatial.distance import pdist, squareform

def create_off_diagonal_matrix(n):
    """
    Creates an n×n matrix where all elements are 0 except for the two diagonals
    adjacent to the main diagonal.
    
    Args:
        n: Size of the square matrix
        
    Returns:
        numpy.ndarray: The resulting matrix
    """
    # Create an n×n matrix of zeros
    matrix = np.zeros((n, n))
    
    # Set the first super-diagonal (above the main diagonal)
    for i in range(n-1):
        matrix[i, i+1] = 1
    
    # Set the first sub-diagonal (below the main diagonal)
    for i in range(n-1):
        matrix[i+1, i] = 1
    
    matrix[0][n-1] = 1
    matrix[n-1][0] = 1
    return matrix

def compute_distance_matrix(positions, alpha):
    """
    Compute the Euclidean distance matrix for a set of positions.
    
    Args:
        positions (list or np.array): A list or array of position tuples/arrays
        
    Returns:
        np.array: A symmetric distance matrix where element [i,j] is the distance
                  between positions[i] and positions[j]
    """
    # Convert to numpy array if not already
    pos_array = np.array(positions)
    n = len(pos_array)
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n, n))
    
    # Compute pairwise distances
    for i in range(n):
        for j in range(n):
            # Euclidean distance between position i and j
            distance_matrix[i, j] = np.linalg.norm(pos_array[i] - pos_array[j])**alpha
    
    return distance_matrix

# Vrf = 500
# V = (12.6, -15.9, 1)
Vrf = 300
V = (38, 36, 21)
NIons = 4
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()
secular_frequencies = [2.0*np.pi*5, 2.0*np.pi*5, 2.0*np.pi*0.1]
equil_positions = [[5e-6, 0, 5e-6],
                   [-5e-6, 0, 5e-6],
                   [-5e-6, 0, -5e-6],
                   [5e-6, 0, -5e-6],
                   ]


Jij_ideal = ig.GlobalInteractionIdeal(NIons)
Jij_ideal = np.array([[0, 1, 1, 0],[0, 0, 0, 1],[0, 0, 0, 1], [0, 0, 0, 0]])
# Jij_ideal = create_off_diagonal_matrix(NIons)
# Jij_ideal = np.array([
#     [0, 1, 1, 1, 0, 0, 0],
#     [1, 0, 0, 1, 1, 0, 0],
#     [1, 0, 0, 1, 0, 1, 0],
#     [1, 1, 1, 0, 1, 1, 1],
#     [0, 1, 0, 1, 0, 0, 1],
#     [0, 0, 1, 1, 0, 0, 1],
#     [0, 0, 0, 1, 1, 1, 0],
# ])
print("secular_frequencies:", secular_frequencies)
print("equil_positions:", equil_positions)



# # Example usage:
# equil_positions = [
#     [-1.85927160e-05, 0.00000000e+00, -6.76215387e-06],
#     [-3.37899564e-05, 0.00000000e+00, 1.79612893e-06],
#     [1.85927160e-05, 0.00000000e+00, -6.76215387e-06],
#     [-1.04893535e-05, 0.00000000e+00, 8.48266971e-06],
#     [0.00000000e+00, 0.00000000e+00, -7.03328954e-06],
#     [1.04893535e-05, 0.00000000e+00, 8.48266971e-06],
#     [3.37899564e-05, 0.00000000e+00, 1.79612893e-06]
# ]

# Jij_ideal = compute_distance_matrix(equil_positions, 1.2)
print("distance_matrix:")
gg.print_matrix(Jij_ideal)

laserset, _ = ig.Jij_Generation(

    Jij_ideal=Jij_ideal, 
    equilibrium_positions=equil_positions, 
    secular_frequencies=secular_frequencies, 
    NIons=NIons
    )

for laser in laserset:
    print(f"frequency:{laser.frequency}, amplitude: {laser.Omegas}")
Jij_exp = gg.Trap_to_Jij(
    equil_positions,
    secular_frequencies,
    NIons,
    laserset
)

print("Jij_ideal:")
gg.print_matrix(Jij_ideal)
print("Jij_exp:")
gg.print_matrix(Jij_exp)
print("laser1:", laserset[0].frequency, laserset[0].Omegas)
inf = ig.Infidelity(Jij_ideal, Jij_exp)
print("infidelity:", inf)
