import numpy as np
import GateGeneration as gg
import Experiment as ex
import MathematicaScripts as scripts
import Circuit as c

Vrf = 500
V = (12.6, -15.9, 1)
NIons = 3
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()

Jij_1 = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
Jij_2 = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]])
c1 = c.Circuit(3)
x_rot = c.Circuit.rotation('x', np.pi)
x_ent = c.Circuit.entanglement('x', Jij_1)
y_ent = c.Circuit.entanglement('y', Jij_2)
c1.addOperator(x_rot)
c1.addOperator(x_ent)
c1.addOperator(y_ent)
c1.addOperator(y_ent)
c1.addOperator(x_rot)
c1.printCircuit()

effective_circuit, fidelity = c1.effective_circuit(secular_frequencies, equil_positions, NIons)
print("fidelity:", fidelity)
c1.print_effective_circuit(effective_circuit)
