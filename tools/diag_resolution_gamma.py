#!/usr/bin/env python3
"""Diagnostic: PL I_max and A_max vs resolution, before/after Lorentzian fix.

For both pyphotonics and qqs Photoluminescence, sweep resolution and gamma
on the diamond test case. Compare pre-fix (`8a5b2b9^`) vs post-fix
(`8a5b2b9`) by git-checking-out the corresponding source files.

Outputs a 2x2 PNG figure (`tools/resolution_gamma_diagnostic.png`):
  row 0: pyphot  | before | after
  row 1: qqs     | before | after
Each panel: I_max and A_max as a function of resolution, one line per gamma.

Run from repo root:
    python tools/diag_resolution_gamma.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path("/home/duguex/pyphotonics")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "qqs" / "lineshape" / "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PRE_COMMIT = "8a5b2b9^"
POST_COMMIT = "8a5b2b9"
RESOLUTIONS = [500, 1000, 4000]
GAMMAS = [0.001, 0.01, 0.05]
DIAMOND = REPO / "test" / "photoluminscence"

CACHE_DIRS = [
    REPO / "__pycache__",
    REPO / "pyphotonics" / "__pycache__",
    REPO / "qqs" / "lineshape" / "src" / "__pycache__",
    REPO / "qqs" / "lineshape" / "src" / "photonics2" / "__pycache__",
]


def clear_cache() -> None:
    for p in CACHE_DIRS:
        if p.exists():
            shutil.rmtree(p, ignore_errors=True)


def checkout(commit: str) -> None:
    subprocess.run(
        ["git", "-C", str(REPO), "checkout", commit, "--",
         "pyphotonics/photoluminescence.py",
         "qqs/lineshape/src/photonics2/photoluminescence.py"],
        check=True, capture_output=True,
    )


def run_pyphot(gamma: float, resolution: int) -> tuple[float, float]:
    from pyphotonics.photoluminescence import Photoluminescence
    p = Photoluminescence(
        ground_state=str(DIAMOND / "CONTCAR_GS"),
        excited_state=str(DIAMOND / "CONTCAR_ES"),
        num_modes=189, method="phonopy",
        phonopy_path=str(DIAMOND / "phonopy"),
        resolution=resolution,
    )
    A, I = p.PL(gamma=gamma, SHR=0, EZPL=0)
    return float(np.max(np.abs(I))), float(np.max(np.abs(A)))


def run_qqs(gamma: float, resolution: int) -> tuple[float, float]:
    from photonics2.photoluminescence import Photoluminescence
    p = Photoluminescence(
        str(DIAMOND / "phonopy" / "band.yaml"),
        "phonopy",
        POSCAR_GRD=str(DIAMOND / "CONTCAR_GS"),
        POSCAR_EX=str(DIAMOND / "CONTCAR_ES"),
        n_defect=1, resolution=resolution,
    )
    p.HuangRhyes()
    p.el_ph(delta_width=6e-3, temperature=0, jtmodes=[])
    p.PL(gamma=gamma, SHR=0, EZPL=0, process="emission")
    p.PLA()
    return float(np.max(np.abs(p.I))), float(np.max(np.abs(p.A)))


def collect(label: str, runner) -> list[dict]:
    rows: list[dict] = []
    for gamma in GAMMAS:
        row: dict = {"gamma": gamma, "I_max": [], "A_max": []}
        for res in RESOLUTIONS:
            i_max, a_max = runner(gamma, res)
            row["I_max"].append(i_max)
            row["A_max"].append(a_max)
            print(f"  {label} gamma={gamma:.3f}  res={res}  I_max={i_max:.4e}  A_max={a_max:.4e}")
        rows.append(row)
    return rows


def main() -> None:
    all_results: dict = {}
    for state, commit in [("before", PRE_COMMIT), ("after", POST_COMMIT)]:
        print(f"\n=== checking out {commit} ===")
        checkout(commit)
        clear_cache()
        print(f"--- {state} pyphot ---")
        all_results[(state, "pyphot")] = collect(state, run_pyphot)
        print(f"--- {state} qqs ---")
        all_results[(state, "qqs")] = collect(state, run_qqs)

    # Always restore post-fix checkout so subsequent runs start fresh.
    checkout(POST_COMMIT)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for row_idx, impl in enumerate(["pyphot", "qqs"]):
        for col_idx, state in enumerate(["before", "after"]):
            ax = axes[row_idx][col_idx]
            rows = all_results[(state, impl)]
            for r in rows:
                ax.semilogy(RESOLUTIONS, r["I_max"], "o-",
                            label=f"I_max γ={r['gamma']:.3f}")
                ax.semilogy(RESOLUTIONS, r["A_max"], "s--",
                            label=f"A_max γ={r['gamma']:.3f}",
                            alpha=0.6)
            ax.set_xlabel("resolution")
            ax.set_ylabel("value")
            ax.set_title(f"{impl} {state}")
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.set_xticks(RESOLUTIONS)

    plt.suptitle("PL I_max and A_max vs resolution — Lorentzian fix comparison")
    plt.tight_layout()
    out = REPO / "tools" / "resolution_gamma_diagnostic.png"
    plt.savefig(out, dpi=120)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()