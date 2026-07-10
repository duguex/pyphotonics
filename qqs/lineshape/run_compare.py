#!/usr/bin/env python3
"""Run test case with both photonics2 versions."""

import os, subprocess, sys

BASE = os.path.dirname(os.path.abspath(__file__))

def run_version(label, script):
    env = os.environ.copy()
    print(f"\n{'='*60}")
    print(f"  Version: {label}")
    print(f"{'='*60}")
    t0 = __import__('time').time()
    result = subprocess.run([sys.executable, "-c", script],
                            cwd=BASE, env=env, capture_output=True, text=True)
    elapsed = __import__('time').time() - t0
    ok = result.returncode == 0
    print(f"  {'✅' if ok else '❌'} Completed in {elapsed:.1f}s")
    for line in result.stdout.splitlines():
        if any(kw in line for kw in ["HuangRhyes", "Error", "OK", "el-ph",
                                      "Start", "skipmode", "Delta_R", "Delta_Q",
                                      "num_modes", "Done", "S_omega", "PL"]):
            print(f"     {line}")
    if result.stderr.strip():
        for line in result.stderr.splitlines()[:5]:
            print(f"  stderr: {line}")
    return ok

# Use Cs3Cu2Br5_STE (smallest band.yaml at 27 MB)
CASE = "Cs3Cu2Br5_STE"
YAML = f"cases/{CASE}/band.yaml"
GRD  = f"cases/{CASE}/GS"
EXC  = f"cases/{CASE}/ES"

# Common preamble to get masses in kg (what the code expects)
preamble = """
from pymatgen.io.vasp.inputs import Poscar
from scipy import constants as const
struct = Poscar.from_file("%s").structure
masses_amu = [s.specie.atomic_mass for s in struct.sites]
masses_kg = [m * const.physical_constants["atomic mass constant"][0] for m in masses_amu]
n_atoms = len(struct)
n_modes = 3 * n_atoms
""" % GRD

script_a = preamble + """
import sys
sys.path.insert(0, 'src')
from photonics2.photoluminescence import Photoluminescence

p = Photoluminescence("%s", "phonopy",
    POSCAR_GRD="%s", POSCAR_EX="%s",
    n_defect=0, resolution=500,
    m=masses_kg, nmodes=n_modes)

hr = p.HuangRhyes()
print(f"HuangRhyes = {hr:.4f}")
print(f"Delta_R = {p.Delta_R:.4f}  Delta_Q = {p.Delta_Q:.4f}")
print(f"skipmodes = {p.skipmodes}  num_modes = {p.numModes}")

p.el_ph(delta_width=6e-3, temperature=0, jtmodes=[])
print(f"S_omega: {len(p.S_omega)} points")
A = p.PL(gamma=1, SHR=0, EZPL=2.0, process="emission")
I = p.PLA()
print(f"PL lineshape: {len(p.I)} points")
print("Done (current)")
""" % (YAML, GRD, EXC)

script_b = preamble + """
import sys
for mod in list(sys.modules):
    if 'photonics2' in mod:
        del sys.modules[mod]
sys.path.insert(0, 'src/lineshape_new_ref')
from photonics2.photoluminescence import Photoluminescence

p = Photoluminescence("%s", "phonopy",
    POSCAR_GRD="%s", POSCAR_EX="%s",
    n_defect=0, resolution=500,
    m=masses_kg, nmodes=n_modes)

hr = p.HuangRhyes()
print(f"HuangRhyes = {hr:.4f}")
print(f"Delta_R = {p.Delta_R:.4f}  Delta_Q = {p.Delta_Q:.4f}")
print(f"skipmodes = {p.skipmodes}  num_modes = {p.numModes}")

p.el_ph(delta_width=6e-3, temperature=0, jtmodes=[])
print(f"S_omega: {len(p.S_omega)} points")
A = p.PL(gamma=1, SHR=0, EZPL=2.0, process="emission")
I = p.PLA()
print(f"PL lineshape: {len(p.I)} points")
print("Done (old vectorized)")
""" % (YAML, GRD, EXC)

ok_a = run_version("Current (src/photonics2/)", script_a)
ok_b = run_version("Old vectorized (src/lineshape_new_ref/)", script_b)

print(f"\n{'='*60}")
print(f"  Summary: A={'✅' if ok_a else '❌'}  B={'✅' if ok_b else '❌'}")
print(f"{'='*60}")
