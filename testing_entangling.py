import MathematicaScripts as scripts
import GateGeneration as gg
import numpy as np

# Set trap parameters
# Vrf = 300
# V = (38, 36, 21)
Vrf = 500
V = (12.6, -15.9, 1)
NIons = 2


#find trap data
secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# secular_frequencies = mathematica_scripts.current_secular_frequencies("ren")
# secular_frequencies = mathematica_scripts.current_secular_frequencies("mid")
print("Secular Frequencies:", secular_frequencies)

equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
# equil_positions = mathematica_scripts.current_equilibrium_positions()
print("Equilibrium Positions:", equil_positions)
# print("ion 1:", equil_positions[0])


mode_frequencies, mode_vectors = gg.eigensystem(equil_positions, secular_frequencies, NIons)
print("mode_frequencies: ", mode_frequencies)
# print("mode_vectors: ", mode_vectors)


laser1 = gg.Laser(frequency=5114170, Omegas=[180e3] * NIons, phase=np.pi/4)
# laser2 = HamiltonianGeneration.Laser(frequency=5108170, Omegas=[180e3] * NIons, phase=np.pi/2)
laserset = [laser1]

total_J_int = gg.Jij(mode_frequencies, mode_vectors, NIons, laserset)

print("Jij matrix: ", total_J_int)

final_Hamiltonian = gg.MS_Hamiltonian(total_J_int, NIons, laserset[0].phase)

print("Final Hamiltonian: ", np.round(final_Hamiltonian, 3))

final_Unitary = gg.unitary_evolution(final_Hamiltonian, .01)

print("Final Unitary:", np.round(final_Unitary, 3))



