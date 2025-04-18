import numpy as np
import Experiment as ex
import GateGeneration as gg
import MathematicaScripts as scripts

I = np.eye(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])
def get_hamiltonian(J, B, c):
    """Constructs the Hamiltonian matrix for the system."""
    H = (-J / 2 * (1 + c) * np.kron(sx, sx)
         - J / 2 * (1 - c) * np.kron(sy, sy)
         - B * np.kron(sz, I)
         - B * np.kron(I, sz))
    return H

def Simulate_H2(J, B, c, t_sim, n):

    # Pauli Matrices
    I = np.eye(2)
    sx = np.array([[0, 1], [1, 0]])
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]])

    # Hamiltonian Construction (H2 problem)
    def get_hamiltonian(J, B, c):
        """Constructs the Hamiltonian matrix for the system."""
        H = (-J / 2 * (1 + c) * np.kron(sx, sx)
             -J / 2 * (1 - c) * np.kron(sy, sy)
             - B * np.kron(sz, I)
             - B * np.kron(I, sz))
        # print("element 1 of hamil:")
        # gg.print_matrix(-2 * J * (1 + c) * np.kron(sx, sx))
        # print("element 1 of unit:")
        # e1 = gg.unitary_evolution(-2 * J * (1 + c) * np.kron(sx, sx), t_sim)
        # gg.print_matrix(e1)
        # print("element 2 of hamil:")
        # gg.print_matrix(- 2 * J * (1 - c) * np.kron(sy, sy))
        # print("element 2 of unit:")
        # e1 = gg.unitary_evolution(- 2 * J * (1 - c) * np.kron(sy, sy), t_sim)
        # gg.print_matrix(e1)    
        print("element 3 of hamil:")
        gg.print_matrix(- B * np.kron(sz, I) - B * np.kron(I, sz))
        print("element 3 of unit:")
        e1 = gg.unitary_evolution(- B * np.kron(sz, I) - B * np.kron(I, sz), t_sim)
        gg.print_matrix(e1)
        return H

    #GOAL HAMILTONIAN
    goal_hamil = get_hamiltonian(J, B, c)
    goal_unitary = gg.unitary_evolution(goal_hamil, t_sim)
    print("Goal Hamiltonian:")
    gg.print_matrix(goal_hamil)
    print("Goal Unitary:")
    gg.print_matrix(goal_unitary)

    Vrf = 500
    V = (12.6, -15.9, 1)
    NIons = 2
    secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
    equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
    # secular_frequencies = scripts.current_secular_frequencies("ren")
    # equil_positions = scripts.current_equilibrium_positions()
    print("secular:", secular_frequencies)
    # print("equil:", equil_positions)
    laser1 = gg.Laser(frequency=5108170, Omegas=[180e3] * NIons, phase=np.pi/2)
    laser2 = gg.Laser(frequency=5108170, Omegas=[180e3] * NIons, phase=0)
    # laser2 = HamiltonianGeneration.Laser(frequency=5108170, Omegas=[180e3] * NIons, phase=np.pi/2)
    laserset_y = [laser1]
    laserset_x = [laser2]


    # print("\n\n\n MS GATE \n\n\n")

    H_MS_x, Jij_x = gg.MS_gate_Hamil(equilibrium_positions=equil_positions, secular_frequencies=secular_frequencies, NIons=NIons, laserset=laserset_x, return_Jij=True)
    print("H_MS_x:")
    gg.print_matrix(H_MS_x, 2)
    H_MS_y, Jij_y = gg.MS_gate_Hamil(equilibrium_positions=equil_positions, secular_frequencies=secular_frequencies, NIons=NIons, laserset=laserset_y, return_Jij=True)
    print("H_MS_y:")
    gg.print_matrix(H_MS_y, 2)

    print("\n\n\n Z GATE \n\n\n")
    z_omega = 5e5
    H_Z = gg.stark_shift([z_omega, z_omega])
    # H_Z = -1 * gg.global_z(2)
    print("H_Z:")
    gg.print_matrix(H_Z)


    #TROTTERIZATION:
    

    t_act_x = -J*(1+c)*t_sim/(H_MS_x[1][2]*n*2)
    t_act_y = -J*(1-c)*t_sim/(H_MS_y[1][2]*n*2)
    t_act_z = -B*2*t_sim/(n*z_omega)
    print("t_act_x:", t_act_x)
    print("t_act_y:", t_act_y)
    print("t_act_z:", t_act_z)


    U_final=np.kron(I, I)
    # Apply Trotter steps
    for _ in range(n):
        U_x = gg.unitary_evolution(H_MS_x, t_act_x)
        U_y = gg.unitary_evolution(H_MS_y, t_act_y)
        U_z = gg.unitary_evolution(H_Z, t_act_z)
        U_final = U_final @ U_x @ U_y @ U_z  # Multiply the unitaries
    print("\n\nTROTTERIZED UNITARIES")
    gg.print_matrix(U_x)
    print()
    gg.print_matrix(U_y)
    print()
    gg.print_matrix(U_z)     
    print("U_final:")
    gg.print_matrix(U_final)
    print("U_goal:")
    gg.print_matrix(goal_unitary)


    diff_matrix = U_final - goal_unitary

    difference = np.linalg.norm(diff_matrix, ord='fro')

    print("difference:", difference)

    return U_final

hbar = 1.054571817e-34
J = 11*hbar
B = 1*hbar
c = 0.2
t_sim = 1

simulating_gate = Simulate_H2(J, B, c, t_sim, 100)

experiment = ex.Experiment(2)
experiment.set_state([0, 1/np.sqrt(2), 1/np.sqrt(2), 0])
experiment.print_state()
experiment.apply_global_gate(simulating_gate)
experiment.print_state()

experiment = ex.Experiment(2)
experiment.set_state([0, -1/np.sqrt(2), 1/np.sqrt(2), 0])
experiment.print_state()
experiment.apply_global_gate(gg.unitary_evolution(get_hamiltonian(J, B, c), t_sim))
experiment.print_state()

print("other eigenvector:")
alpha = np.sqrt(4*(B**2) + (J**2)*(c**2))
stateEntry1 = np.sqrt((alpha+2*B)/(2*alpha))
stateEntry2 = np.sqrt((alpha-2*B)/(2*alpha))
print("stateEntry1:", stateEntry1)
print("stateEntry2:", stateEntry2)

experiment = ex.Experiment(2)
experiment.set_state([stateEntry1, 0, 0, stateEntry2])
experiment.print_state()
experiment.apply_global_gate(simulating_gate)
experiment.print_state()

experiment = ex.Experiment(2)
experiment.set_state([stateEntry1, 0, 0, stateEntry2])
experiment.print_state()
experiment.apply_global_gate(gg.unitary_evolution(get_hamiltonian(J, B, c), t_sim))
experiment.print_state()
