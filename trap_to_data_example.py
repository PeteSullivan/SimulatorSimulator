import MathematicaScripts as scripts
import GateGeneration as gg

# Set trap parameters
Vrf = 300
V = (38, 36, 21)
NIons = 5

#find secular frequencies
secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# secular_frequencies = mathematica_scripts.current_secular_frequencies("ren")
# secular_frequencies = mathematica_scripts.current_secular_frequencies("mid")
print("Secular Frequencies:", secular_frequencies)


#find equilibrium positions
equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
# equil_positions = mathematica_scripts.current_equilibrium_positions()
# print("Equilibrium Positions:", equil_positions)

#find mode frequencies and vectors
mode_frequencies, mode_vectors = gg.eigensystem(equil_positions, secular_frequencies, NIons)
print("mode_frequencies: ", mode_frequencies)
print("mode_vectors: ", mode_vectors)
