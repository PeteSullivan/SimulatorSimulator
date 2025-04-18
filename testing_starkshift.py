import GateGeneration as gg

Omegas = [6, 6]

Z_gate = gg.stark_shift(Omegas=Omegas)

print("Z_gate:")
gg.print_matrix(Z_gate)