import MathematicaScripts as scripts
import numpy as np
import GateGeneration as gg

Vrf = 500
V = (12.6, -15.9, 1)
NIons = 4
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()
print("secular:", secular_frequencies)
print("equil:", equil_positions)
_, NormalModeEigVecs = gg.eigensystem(equil_positions, secular_frequencies, NIons)
JDes = np.array([
    [0, 1, 1, 1],
    [1, 0, 1, 1],
    [1, 1, 0, 1],
    [1, 1, 1, 1]
    ])

weights = scripts.Jij_to_weights(JDes, NormalModeEigVecs)

# gg.print_matrix(weights)
print("weights:", weights)