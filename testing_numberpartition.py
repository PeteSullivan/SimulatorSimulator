import numpy as np
import Experiment as ex
import GateGeneration as gg
import MathematicaScripts as scripts
import InteractionGeneration as ig

def ZZ_interaction(set):
    NIons = len(set)
    ZZ = np.zeros([NIons, NIons], dtype="complex")
    for i in range(NIons):
        for j in range(NIons):
            if i != j:
                ZZ[i][j] = 2 * set[i] * set[j]
    return ZZ

set = [1, 2, 2, 3]
NIons = len(set)

ZZ = ZZ_interaction(set)
gg.print_matrix(ZZ)

trotterization_steps = 10

experiment = ex.Experiment(NIons)

experiment.set_state([1/4]*16)
experiment.print_state()
experiment.global_entanglement(20, ZZ, 'z', NIons)
# for i in range(trotterization_steps*10):
#     experiment.global_entanglement(np.pi/trotterization_steps, ZZ, 'z', NIons)
#     # experiment.global_rotation(np.pi/trotterization_steps, 'x', NIons)
#     # experiment.print_state()
experiment.print_state()
experiment.global_rotation(0.6, 'x', NIons)
experiment.print_state()
experiment.measure()