"""Example: Shifted NV center in diamond."""
from pyphotonics import Photoluminescence
import numpy as np
import matplotlib.pyplot as plt

masses = np.zeros(63)
for i in range(63):
    masses[i] = 12.011 * 1.660539040e-27
masses[62] = 14.007 * 1.660539040e-27

p = Photoluminescence(
    ground_state="MOD_CONTCAR_GS",
    excited_state="MOD_CONTCAR_ES",
    num_modes=189,
    method="phonopy",
    phonopy_path="phonopy/",
    masses=masses,
)

print("Delta_R=", p.Delta_R)
print("Delta_Q=", p.Delta_Q)
print("HuangRhys=", p.HuangRhys)
print("DebyeWaller=", np.exp(-p.HuangRhys))

plt.figure(figsize=(10, 10))
plt.plot(p.S_omega)
plt.ylabel("$S(\hbar\omega)$")
plt.xlabel("Phonon energy (meV)")
plt.xlim(0, 200)
plt.savefig("S_omega", bbox_inches="tight")

p.write_S("S")

A, I = p.PL(2, 2, 1.95)

plt.figure(figsize=(10, 10))
plt.plot(np.abs(I))
plt.ylabel("$I(\hbar\omega)$")
plt.xlabel("Photon energy (eV)")
plt.xlim(1200, 2000)
x_values, labels = plt.xticks()
labels = [float(x) / p.resolution for x in x_values]
plt.xticks(x_values, labels)
plt.ylim(0, 600)
plt.savefig("I", bbox_inches="tight")
