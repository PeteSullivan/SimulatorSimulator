import GateGeneration as gg
import Experiment as ex
import numpy as np
import InteractionGeneration as intg
class Circuit:

    class rotation:
        def __init__(self, direction, angle, t):
            self.angle = angle
            self.direction = direction


    class entanglement:
        def __init__(self, direction, Jij, t):
            self.direction = direction
            self.Jij = Jij


    def __init__(self, NIons):
        NIons = NIons
        self.operators = []

    def addOperator(self, operator):
        self.operators.append(operator)

    def printCircuit(self):
        for operator in self.operators:

            if type(operator) == self.rotation:
                print(f"rotation around {operator.direction} by {operator.angle/np.pi} pi")
            elif type(operator) == self.entanglement:
                print(f"entanglement around {operator.direction}, interactions:")
                gg.print_matrix(operator.Jij)

    def effective_circuit(self, secular_frequencies, equil_positions, NIons):
        laser_params = []
        total_fidelity = 1
        time = 0
        for operator in self.operators:
            if type(operator) == self.rotation:
                if operator.direction == 'x':
                    runtime = operator.angle / 180e3
                    laser = gg.Laser(frequency=12.6e9, Omegas=[180e3]*NIons, phase=0)
                    laser_param = [laser, [time, time + runtime]]
                    laser_params.append(laser_param)
                    time += runtime + 1e-9

                elif operator.direction == 'y':
                    runtime = operator.angle / 180e3
                    laser = gg.Laser(frequency=12.6e9, Omegas=[180e3]*NIons, phase=np.pi/2)
                    laser_param = [laser, [time, time + runtime]]
                    laser_params.append(laser_param)
                    time += runtime + 1e-9

                else:
                    runtime = operator.angle / 180e3
                    laser = gg.Laser(frequency=1, Omegas=1)
                    laser_param = [laser, time, time + runtime]
                    laser_params.append(laser_param)
                    time += runtime + 1e-9

            else:
                laserset, infidelity = intg.Jij_Generation(
                    equilibrium_positions=equil_positions,
                    secular_frequencies=secular_frequencies,
                    NIons=NIons,
                    Jij_ideal=operator.Jij
                )
                for laser in laserset:
                    runtime = 1e-6
                    laser_param = [laser, [time, time + runtime]]
                    laser_params.append(laser_param)
                time += runtime + 1e-9
                total_fidelity = total_fidelity * (1 - infidelity)

        return laser_params, total_fidelity

    def print_effective_circuit(self, effective_circuit):
        for laser_param in effective_circuit:
            laser = laser_param[0]
            print(f"laser with frequency {laser.frequency} from {np.round(laser_param[1][0], 9)} to {np.round(laser_param[1][1], 9)}")