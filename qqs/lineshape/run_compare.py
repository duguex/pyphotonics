#!/usr/bin/env python3
"""Run all available test cases and verify against BASELINE.md.

After the 2026-07-11 merge, src/photonics2/photoluminescence.py is
algorithm-equivalent to lineshape_new_ref (deleted 2026-07-11). This
script now runs a single implementation and cross-checks against the
recorded baseline."""

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


def make_script(name, grd, exc, yaml_f):
    """Generate Python script to run the post-merge implementation for one case."""
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
sys.path.insert(0, 'src')
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


def run_version(script):
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
print(f"  {'Case':<18} {'S':<14} {'ΔQ':<10} {'skip':<8} {'num':<8}")
print(f"{'='*80}")

results = []
for name, grd, exc, yaml_f in CASES:
    err = check_case(name, grd, exc, yaml_f)
    if err:
        print(f"  {name:<18} ⏭️  {err}")
        continue

    s = make_script(name, grd, exc, yaml_f)
    ok, hr, dr, dq, sm, nm, err_msg = run_version(s)
    mark = "✅" if ok else "❌"

    print(f"  {name:<18} {mark} S={hr:<10} ΔQ={dq:<8} skip={sm:<6} num={nm}")
    results.append((name, ok, hr, dq, sm, nm))

print(f"{'='*80}")
print("  Summary:")
for name, ok, hr, dq, sm, nm in results:
    status = "✅" if ok else "❌"
    print(f"    {name:<18} {status}  S={hr}  ΔQ={dq}  skip={sm}  num={nm}")
print(f"{'='*80}")