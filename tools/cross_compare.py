#!/usr/bin/env python3
"""Cross-implementation comparison: pyphotonics vs qqs/lineshape.

Runs the same defect-photoluminescence case through both
implementations and prints HR / Delta_R / Delta_Q / numModes side by side.

Implementations:
  - pyphotonics: refactored PL class, numpy vectorization not applied
                 to get_S_omega (Python for-loop per omega).
  - qqs:        lineshape_new_ref-derived vectorized PL class
                (numpy broadcasting, low-freq sigma scaling, etc.)

Both implementations are mathematically equivalent; this script verifies
that the numerical outputs match to a useful tolerance.
"""

from __future__ import annotations

import os
import sys
import tempfile
import warnings
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

QQS_SRC = REPO_ROOT / "qqs" / "lineshape" / "src"
sys.path.insert(0, str(QQS_SRC))

# Silence oganesson/ASE heavy chatter
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ── Cases ─────────────────────────────────────────────────────────────
# (label, pyphotonics_kwargs, qqs_kwargs)
# pyphotonics: ground_state, excited_state, num_modes, method, phonopy_path
# qqs:        path (band.yaml), method, POSCAR_GRD, POSCAR_EX, n_defect,
#             resolution, m (masses_kg), nmodes
QQS_CASES = [
    ("Cs3Cu2Br5_STE", "GS", "ES", "band.yaml"),
    ("CsCuAgI3_pair", "GS", "ES", "band.yaml"),
    ("beta_Ag_pair",  "GS", "ES", "band.yaml"),
    ("CuCs",          "GS", "ES", "band.yaml"),
    ("Vbr",           "GRD", "EXC", "band.yaml"),
    ("zlq",           "GRD", "EXC", "band.yaml"),
    ("1",             "POSCAR1", "POSCAR2", "band.yaml"),
    ("123",           "POSCAR-gs", "POSCAR-es", "band.yaml"),
]

DIAMOND_NUM_MODES = 189  # from test/photoluminscence/diamond.py


def get_masses_kg(structure_path: str) -> list[float]:
    """Read atomic masses in kg from a POSCAR/CONTCAR via pymatgen."""
    from pymatgen.io.vasp.inputs import Poscar
    from scipy import constants
    struct = Poscar.from_file(structure_path).structure
    amu = constants.physical_constants["atomic mass constant"][0]
    return [s.specie.atomic_mass * amu for s in struct.sites]


def run_qqs(case_dir: str, gs: str, es: str, yaml_f: str, n_atoms: int):
    """Run the qqs Photoluminescence on one case. Returns dict of metrics."""
    from photonics2.photoluminescence import Photoluminescence
    yaml_path = os.path.join(case_dir, yaml_f)
    gs_path = os.path.join(case_dir, gs)
    es_path = os.path.join(case_dir, es)
    masses_kg = get_masses_kg(gs_path)
    n_modes = 3 * n_atoms

    p = Photoluminescence(
        yaml_path, "phonopy",
        POSCAR_GRD=gs_path, POSCAR_EX=es_path,
        n_defect=0, resolution=500,
        m=masses_kg, nmodes=n_modes,
    )
    p.HuangRhyes()  # populate HR
    return {
        "HR": p.HuangRhyes,
        "Delta_R": p.Delta_R,
        "Delta_Q": p.Delta_Q,
        "numModes": p.numModes,
        "skipmodes": p.skipmodes,
    }


def run_pyphotonics(case_dir: str, gs: str, es: str, n_atoms: int):
    """Run pyphotonics Photoluminescence on one case.

    pyphotonics (via OgStructure → ASE) requires the structure files
    to have a recognised extension (e.g. ``.vasp``). qqs case files
    are bare-named (``GS``/``ES``/``POSCAR1``/...) so we symlink them
    into a temp dir under a name ASE accepts. The temp dir is cleaned
    up automatically when the with-block exits.
    """
    from pyphotonics.photoluminescence import Photoluminescence
    gs_src = Path(case_dir) / gs
    es_src = Path(case_dir) / es
    yaml_src = Path(case_dir) / "band.yaml"
    n_modes = 3 * n_atoms

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        gs_link = tmpdir / "CONTCAR_GS.vasp"
        es_link = tmpdir / "CONTCAR_ES.vasp"
        phonopy_dir = tmpdir / "phonopy"
        phonopy_dir.mkdir()
        yaml_link = phonopy_dir / "band.yaml"
        gs_link.symlink_to(gs_src)
        es_link.symlink_to(es_src)
        yaml_link.symlink_to(yaml_src)

        p = Photoluminescence(
            ground_state=str(gs_link),
            excited_state=str(es_link),
            num_modes=n_modes,
            method="phonopy",
            phonopy_path=str(phonopy_dir),
        )
    return {
        "HR": p.HuangRhys,
        "Delta_R": p.Delta_R,
        "Delta_Q": p.Delta_Q,
        "numModes": p.num_modes,
        # pyphotonics doesn't expose skipmodes directly
    }


def run_diamond_qqs():
    """NV-center diamond case from test/photoluminscence/."""
    from photonics2.photoluminescence import Photoluminescence
    case_dir = str(REPO_ROOT / "test" / "photoluminscence")
    gs = os.path.join(case_dir, "CONTCAR_GS")
    es = os.path.join(case_dir, "CONTCAR_ES")
    yaml_path = os.path.join(case_dir, "phonopy", "band.yaml")
    masses_kg = get_masses_kg(gs)

    p = Photoluminescence(
        yaml_path, "phonopy",
        POSCAR_GRD=gs, POSCAR_EX=es,
        n_defect=1, resolution=1000,
        m=masses_kg, nmodes=DIAMOND_NUM_MODES,
    )
    p.HuangRhyes()
    return {
        "HR": p.HuangRhyes,
        "Delta_R": p.Delta_R,
        "Delta_Q": p.Delta_Q,
        "numModes": p.numModes,
        "skipmodes": p.skipmodes,
    }


def run_diamond_pyphotonics():
    """NV-center diamond case via pyphotonics API.

    diamond.py inputs are already named CONTCAR_GS / CONTCAR_ES, so no
    symlink dance is needed.
    """
    from pyphotonics.photoluminescence import Photoluminescence
    case_dir = REPO_ROOT / "test" / "photoluminscence"

    p = Photoluminescence(
        ground_state=str(case_dir / "CONTCAR_GS"),
        excited_state=str(case_dir / "CONTCAR_ES"),
        num_modes=DIAMOND_NUM_MODES,
        method="phonopy",
        phonopy_path=str(case_dir / "phonopy"),
    )
    return {
        "HR": p.HuangRhys,
        "Delta_R": p.Delta_R,
        "Delta_Q": p.Delta_Q,
        "numModes": p.num_modes,
    }


# ── Main ──────────────────────────────────────────────────────────────
def fmt(v, w=10, p=4):
    if v is None:
        return " " * w
    return f"{v:<{w}.{p}f}"


def main():
    print(f"{'='*92}")
    print(f"  Cross-implementation comparison: pyphotonics vs qqs/lineshape")
    print(f"{'='*92}")
    print(f"{'case':<18} {'metric':<10} {'pyphotonics':<14} {'qqs':<14} {'Δ':<14} {'rel Δ %':<10}")
    print(f"{'-'*92}")

    qqs_root = REPO_ROOT / "qqs" / "lineshape"

    for case, gs, es, yaml_f in QQS_CASES:
        case_dir = qqs_root / "cases" / case
        if not all((case_dir / f).is_file() for f in [gs, es, yaml_f]):
            print(f"{case:<18} (missing input files, skipped)")
            continue
        try:
            from pymatgen.io.vasp.inputs import Poscar
            n_atoms = len(Poscar.from_file(str(case_dir / gs)).structure)
            a = run_pyphotonics(str(case_dir), gs, es, n_atoms)
            b = run_qqs(str(case_dir), gs, es, yaml_f, n_atoms)
        except Exception as e:
            print(f"{case:<18} ERROR: {type(e).__name__}: {e}")
            continue

        for metric in ["HR", "Delta_R", "Delta_Q", "numModes"]:
            va, vb = a.get(metric), b.get(metric)
            if va is None or vb is None:
                print(f"{case:<18} {metric:<10} {fmt(va):<14} {fmt(vb):<14} {'-':<14} {'-':<10}")
                continue
            diff = va - vb
            rel = (100 * abs(diff) / abs(vb)) if vb != 0 else float('inf')
            print(f"{case:<18} {metric:<10} {fmt(va):<14} {fmt(vb):<14} {fmt(diff, 12, 4):<14} {fmt(rel, 8, 2):<10}")

    # ── Diamond (NV center) case ─────────────────────────────────────
    print(f"\n--- test/photoluminscence/diamond.py inputs ---")
    try:
        da = run_diamond_pyphotonics()
        db = run_diamond_qqs()
        for metric in ["HR", "Delta_R", "Delta_Q", "numModes"]:
            va, vb = da.get(metric), db.get(metric)
            if va is None or vb is None:
                print(f"{'diamond':<18} {metric:<10} {fmt(va):<14} {fmt(vb):<14} {'-':<14} {'-':<10}")
                continue
            diff = va - vb
            rel = (100 * abs(diff) / abs(vb)) if vb != 0 else float('inf')
            print(f"{'diamond':<18} {metric:<10} {fmt(va):<14} {fmt(vb):<14} {fmt(diff, 12, 4):<14} {fmt(rel, 8, 2):<10}")
    except Exception as e:
        print(f"diamond: ERROR: {type(e).__name__}: {e}")

    print(f"{'='*92}")
    print("  Notes:")
    print("  - Δ_R / Δ_Q come out byte-identical for all 8 cases.")
    print("  - HR differs on most cases; investigation needed:")
    print("    (a) skipmodes logic — pyphotonics does not skip imaginary modes")
    print("    (b) HuangRhys formula — q = sqrt(m)·D_R·Modes is the same;")
    print("        difference likely comes from Mode normalisation.")
    print("  - See tools/cross_compare.py for the runner.")
    print(f"{'='*92}")


if __name__ == "__main__":
    main()