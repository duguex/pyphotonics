#!/usr/bin/env python3
"""Run all available test cases with both photonics2 versions."""

import os, subprocess, sys

BASE = os.path.dirname(os.path.abspath(__file__))

# Case definitions: (name, grd_file, exc_file, yaml_file)
# grd/exc are relative to cases/<name>/
CASES = [
    ("Cs3Cu2Br5_STE", "GS", "ES", "band.yaml"),
    ("CsCuAgI3_pair", "GS", "ES", "band.yaml"),
    ("beta_Ag_pair",   "GS", "ES", "band.yaml"),
    ("CuCs",           "GS", "ES", "band.yaml"),
    ("Vbr",            "GRD", "EXC", "band.yaml"),
    ("zlq",            "GRD", "EXC", "band.yaml"),
    ("1",              "POSCAR1", "POSCAR2", "band.yaml"),
    ("123",            "POSCAR-gs", "POSCAR-es", "band.yaml"),
]


def check_case(name, grd, exc, yaml_f):
    """Verify files exist on disk."""
    d = os.path.join(BASE, "cases", name)
    ok = all(os.path.isfile(os.path.join(d, f)) for f in [grd, exc, yaml_f])
    if not ok:
        missing = [f for f in [grd, exc, yaml_f] if not os.path.isfile(os.path.join(d, f))]
        return f"missing: {missing}"
    return None


def make_script(name, grd, exc, yaml_f, version_path):
    """Generate Python script to run one version for one case."""
    yaml = f"cases/{name}/{yaml_f}"
    g = f"cases/{name}/{grd}"
    e = f"cases/{name}/{exc}"
    return f"""
from pymatgen.io.vasp.inputs import Poscar
from scipy import constants as const
struct = Poscar.from_file("{g}").structure
masses_amu = [s.specie.atomic_mass for s in struct.sites]
masses_kg = [m * const.physical_constants["atomic mass constant"][0] for m in masses_amu]
n_atoms = len(struct)
n_modes = 3 * n_atoms

import sys
sys.path.insert(0, '{version_path}')
from photonics2.photoluminescence import Photoluminescence

p = Photoluminescence("{yaml}", "phonopy",
    POSCAR_GRD="{g}", POSCAR_EX="{e}",
    n_defect=0, resolution=500,
    m=masses_kg, nmodes=n_modes)

hr = p.HuangRhyes()
print(f"HuangRhyes={{hr:.4f}}")
print(f"Delta_R={{p.Delta_R:.4f}}  Delta_Q={{p.Delta_Q:.4f}}")
print(f"skipmodes={{p.skipmodes}}  num_modes={{p.numModes}}")

p.el_ph(delta_width=6e-3, temperature=0, jtmodes=[])
print(f"S_omega:{{len(p.S_omega)}} points")
A = p.PL(gamma=1, SHR=0, EZPL=2.0, process="emission")
I = p.PLA()
print(f"PL:{{len(p.I)}} points")
print("OK")
"""

def run_version(label, script):
    env = os.environ.copy()
    try:
        result = subprocess.run([sys.executable, "-c", script],
                                cwd=BASE, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return False, "TIMEOUT", "", "", "", "", ""
    ok = result.returncode == 0
    hr = dr = dq = sm = nm = "?"
    for line in result.stdout.splitlines():
        if "HuangRhyes=" in line:
            hr = line.split("=")[-1].strip()
        elif "Delta_R" in line:
            for tok in line.split():
                if tok.startswith("Delta_R="): dr = tok.split("=")[1]
                if tok.startswith("Delta_Q="): dq = tok.split("=")[1]
        elif "skipmodes" in line:
            for tok in line.split():
                if tok.startswith("skipmodes="): sm = tok.split("=")[1]
                if tok.startswith("num_modes="): nm = tok.split("=")[1]
    err = result.stderr.strip()[:100] if result.stderr else ""
    return ok, hr, dr, dq, sm, nm, err


print(f"{'='*80}")
print(f"  {'Case':<18} {'A(HR)':<14} {'B(HR)':<14} {'A(skip)':<10} {'B(skip)':<10} {'A(ΔQ)':<10} {'B(ΔQ)':<10}")
print(f"{'='*80}")

results = []
for name, grd, exc, yaml_f in CASES:
    err = check_case(name, grd, exc, yaml_f)
    if err:
        print(f"  {name:<18} ⏭️  {err}")
        continue

    sa = make_script(name, grd, exc, yaml_f, "src")
    sb = make_script(name, grd, exc, yaml_f, "src/lineshape_new_ref")

    ok_a, hr_a, dr_a, dq_a, sm_a, nm_a, err_a = run_version("A", sa)
    ok_b, hr_b, dr_b, dq_b, sm_b, nm_b, err_b = run_version("B", sb)

    a_mark = "✅" if ok_a else "❌"
    b_mark = "✅" if ok_b else "❌"
    hr_match = "✓" if ok_a and ok_b and hr_a == hr_b else "✗"

    print(f"  {name:<18} {a_mark} S={hr_a:<10} {b_mark} S={hr_b:<10} "
          f"skip={sm_a:<6} skip={sm_b:<6} ΔQ={dq_a:<6} ΔQ={dq_b:<6} {hr_match}")

    results.append((name, ok_a, ok_b, hr_a, hr_b, hr_a == hr_b if ok_a and ok_b else False))

print(f"{'='*80}")
print("  Summary:")
for name, ok_a, ok_b, hr_a, hr_b, match in results:
    status = "✅" if ok_a and ok_b else "❌"
    match_s = " ✓" if match else " ✗ MISMATCH"
    print(f"    {name:<18} {status}  A:{hr_a}  B:{hr_b}{match_s}")
print(f"{'='*80}")
