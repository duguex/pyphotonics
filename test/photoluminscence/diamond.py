"""Example: Photoluminescence calculation for NV defect in diamond."""
from pyphotonics import Photoluminescence
import numpy as np
import matplotlib.pyplot as plt

phonopy_path = "phonopy/"
p = Photoluminescence(
    ground_state="CONTCAR_GS",
    excited_state="CONTCAR_ES",
    num_modes=189,
    method="phonopy",
    phonopy_path=phonopy_path,
)

print("Delta_R=", p.Delta_R)
print("Delta_Q=", p.Delta_Q)
print("HuangRhys=", p.HuangRhys)

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
