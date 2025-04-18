import MathematicaScripts as scripts
import Circuit as c
import numpy as np
import GateGeneration as gg


#get trap information
Vrf = 500
V = (12.6, -15.9, 1)
NIons = 4
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()
print("secular:", secular_frequencies)
print("equil:", equil_positions)
mode_frequencies, mode_vectors = gg.eigensystem(equil_positions, secular_frequencies, NIons)

#define a JDes to simulate
JDes = np.array([
    [0, 1, 1, 1],
    [1, 0, 1, 1],
    [1, 1, 0, 1],
    [1, 1, 1, 1]
    ])

#find weights to simulate JDes
weights = scripts.Jij_to_weights(JDes, mode_vectors)

#build circuit to simulate JDes
# gg.print_matrix(weights)
print("weights:", weights)

cir = c.Circuit(NIons)

z_rot = c.Circuit.rotation('z', np.pi)

Jij_1 = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
Jij_2 = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]])
c1 = c.Circuit(3)
x_rot = c.Circuit.rotation('x', np.pi, 0.05)
x_ent = c.Circuit.entanglement('x', Jij_1, 0.05)
y_ent = c.Circuit.entanglement('y', Jij_2, 0.05)
c1.addOperator(x_rot)
c1.addOperator(x_ent)
c1.addOperator(y_ent)
c1.addOperator(y_ent)
c1.addOperator(x_rot)
c1.printCircuit()

effective_circuit, fidelity = c1.effective_circuit(secular_frequencies, equil_positions, NIons)
print("fidelity:", fidelity)
c1.print_effective_circuit(effective_circuit)