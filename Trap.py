import MathematicaScripts as scripts


class Trap():
    def __init__(self, name, NIons, V, Vrf, currentData=False):
        self.name = name
        self.NIons = NIons
        if currentData:
            self.secular_frequencies = scripts.current_secular_frequencies(name)
            self.equil_positions = scripts.current_equilibrium_positions()
        else:
            self.secular_frequencies = scripts.find_secular_frequencies("ren", NIons, Vrf, V)
            self.equil_positions = scripts.find_equilibrium_positions(NIons, self.secular_frequencies)
