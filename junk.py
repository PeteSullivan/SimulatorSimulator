import GateGeneration as gg
import numpy as np
import Experiment as ex

Jij = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
Jij2 = np.array([[1, 1, 1],[1,1,1],[1,1,1]])
NIons = 20
e = ex.Experiment(NIons)
e.global_rotation(np.pi/2, 'y', NIons)
e.print_state()
e.global_entanglement('z', Jij, np.pi/4, NIons)
e.print_state()
e.global_entanglement('x', Jij, np.pi/4, NIons)
e.print_state()
e.global_rotation(-np.pi/2, 'y', NIons)
e.print_state()
