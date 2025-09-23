import MathematicaScripts as scripts
import Circuit as c
import numpy as np
import GateGeneration as gg
import InteractionGeneration as ig

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
# mode_vectors = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, -1, -1, 1], [1, -1, 1, -1]])
print("mode vectors:", mode_vectors)
#define a JDes to simulate

pi4 = np.pi / 4
JDes = np.array([
    [0, 0.5, 0.5, 0],
    [0, 0, 0, 0.5],
    [0, 0, 0, 0.5],
    [0, 0, 0, 0]
    ])


laserset, inf = ig.Jij_Generation_NEW(equil_positions, secular_frequencies, 4, JDes)
print("inf:", inf)

# #find weights to simulate JDes
# weights = scripts.Jij_to_weights(JDes, mode_vectors)

# #build circuit to simulate JDes
# # gg.print_matrix(weights)
# print("weights:", weights)

# Jij = gg.Jij_from_weights(mode_frequencies=mode_frequencies, mode_vectors=mode_vectors, NIons=NIons, weights=weights)
# print("Generated Jij:")
# gg.print_matrix(Jij)
# cir = c.Circuit(NIons)
