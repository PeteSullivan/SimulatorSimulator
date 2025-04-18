import GateGeneration as gg
import Experiment as Exp
import MathematicaScripts as scripts
import numpy as np

#find trap data
Vrf = 500
V = (12.6, -15.9, 1)
NIons = 2
# secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
# equil_positions = scripts.find_equilibrium_positions(NIons, secular_frequencies)
secular_frequencies = scripts.current_secular_frequencies("ren")
equil_positions = scripts.current_equilibrium_positions()
print("Secular Frequencies:", secular_frequencies)

#initialize experiment
experiment = Exp.Experiment(2)
experiment.print_state()


laser1 = gg.Laser(frequency=5114170, Omegas=[12.6e9] * NIons, phase=np.pi/2)
# laser2 = HamiltonianGeneration.Laser(frequency=5108170, Omegas=[180e3] * NIons, phase=np.pi/2)
laserset = [laser1]
time = 4e-6
gate1 = gg.MS_gate(equilibrium_positions=equil_positions, secular_frequencies=secular_frequencies, NIons=NIons, laserset=laserset, time=480e-6)

# print("MS gate 1:", gate1)

experiment.apply_global_gate(gate1)
experiment.print_state()


"""
progress: unitaries seem good, sqr is working well, experiment is generally working
to do:
test MS gate with accurate omegas, make sure Experiment is accurate not rounding or anything, 
see how multiple gates can simulate the H2 problem, bring results/questions for tuesday

long term- given a Jij, generate the laserset to simulate it (phil's code for global multimodal?)
given a hamiltonian, can it be simulated?
new laser class, initialize as global, precise, or gaussian, return different omega dependencies
refine sqr class for different states (new frequency goal thing), z-rotations (Ry(pi/2)Rx(pi)Ry(-pi/2) maybe)
format experiments as X-rotations, Y-rotations, Z-rotations, entangling, measurements, repeat
"""
