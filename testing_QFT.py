import GateGeneration as gg
import numpy as np
import MathematicaScripts as scripts

Vrf = 500
V = (12.6, -15.9, 1)
NIons = 2
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()
print("Secular Frequencies:", secular_frequencies)

