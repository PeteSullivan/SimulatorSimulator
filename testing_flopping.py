import GateGeneration as gg
import numpy as np
NIons = 1
piFlip = np.pi / 200e3 #pi / Omega for full rotation
z_rotation = 1
h=6.62607015e-34

laser = gg.Laser(frequency=12.6e9, Omegas=[200e3], phase=0)
H1 = gg.sqr(t=piFlip, laser=laser)
print("phase:", laser.phase, "\nH:", H1)
unitary1 = gg.unitary_evolution(H=H1, t=0.5)
print("Unitary:", unitary1)

laser = gg.Laser(frequency=12.6e9, Omegas=[200e3], phase=np.pi/2)
H2 = gg.sqr(t=piFlip, laser=laser)
print("phase:", laser.phase, "\nH:", H2)
unitary2 = gg.unitary_evolution(H=H2, t=piFlip)
print("Unitary:", unitary2)

laser = gg.Laser(frequency=12.6e9, Omegas=[200e3], phase=0)
H3 = gg.sqr(t=piFlip, laser=laser)
print("phase:", laser.phase, "\nH:", H3)
unitary3 = gg.unitary_evolution(H=H3, t=7*piFlip/2)
print("Unitary:", unitary3)

H_final = unitary1 @ unitary2 @ unitary3
print("full rotation:", np.round(H_final, 4))
