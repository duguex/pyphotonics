# qqs/lineshape/src/photonics2/ — Dormant Modules Development Notes

**Date**: 2026-07-11
**Context**: After the 2026-07-11 photoluminescence merge
(`docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`),
the active PL code in `photoluminescence.py` / `plott.py` is the
numpy-vectorized implementation from `lineshape_new_ref/`. This directory
additionally contains 10 modules that are not part of the active PL
pipeline. They are **dormant** — kept on disk for potential future
development direction (Jahn-Teller physics, defect embedding, harmonic
analysis).

## Active vs. dormant

| File | Status | Role |
|------|--------|------|
| `photoluminescence.py` | **active** | Core PL computation (Huang-Rhys, S(ω), PL line-shape) |
| `plott.py` | **active** | Plotting helpers (post-merge uses lineshape_new_ref version) |
| `configuration_coordinate.py` | **active** | POSCAR reading helper used by `photoluminescence.py` |
| `xyz.py` | **active** | XYZ-format reading helper |
| `__init__.py` | n/a | Empty package marker |
| `jt.py`, `jtsoc.py`, `ht.py`, `schrodinger.py`, `embedding.py`, `ebedding.py`, `hermite.py`, `nu.py`, `un.py`, `constants.py` | **dormant** | See below |

## Dormant modules

All ten dormant modules have a `DEPRECATED/DORMANT` docstring banner at
the top of the file describing their status. The 2026-07-11 cleanup
applied the following fixes:

- **`jt.py`** (773 lines): import path `photonics.ht` → `photonics2.ht`,
  `photonics.hermite` → `photonics2.hermite`. Now imports cleanly.
- **`jtsoc.py`** (1265 lines): import path `photonics.xyz` →
  `photonics2.xyz`. Now imports cleanly.
- **`embedding.py`** (171 lines): removed extra quote in
  `"defect_force_constant"",0)` (line 134); added `pass` to empty
  `else:` branch (line 152). Now imports cleanly.
- **`nu.py`** (11 lines): removed module-level `print()` statements that
  fired on import. Now imports cleanly.
- **`un.py`** (3 lines): removed module-level `print()` statements.
  File now contains only an `import numpy as np`.
- **`ebedding.py`**, **`ht.py`**, **`schrodinger.py`**, **`hermite.py`**,
  **`constants.py`**: untouched, were already import-clean.

### Module catalogue

| Module | Lines | Planned role | Notes |
|--------|------:|--------------|-------|
| `jt.py` | 773 | Jahn-Teller (E⊗e) solver — chiral phonon basis, original-basis split E, T↔E conversion, `MnJT`, `Spectra` | dormant; depends on `sympy` / `joblib` (installed in this env); no callers |
| `jtsoc.py` | 1265 | JT + spin-orbit coupling (SOC) — adds `_soc` variants of `TtosplitE_mp`, `splitEtoT_mp`, `onlyEmode` | dormant; `aps` class with `Cs(T)`, `temp_soc`; no callers |
| `ht.py` | 284 | Huang-Rhys / harmonic analysis (class `phonon`, `band`, fns `htt`, `odd`) | dormant; no callers |
| `schrodinger.py` | 61 | 1D Schrödigner solver class `Schrodinger` | dormant; **functionality overlaps with `carriercapture/potential.py`** which is the active Schrödinger solver for the main `pyphotonics` package. If JT work needs 1D Schrödigner, prefer extending `carriercapture/potential.py` over reviving this. |
| `embedding.py` | 171 | Defect-host lattice embedding — `embedding` class, `conference_matrix` | dormant; pymatgen-based; no callers |
| `ebedding.py` | 7 | Empty placeholder (typo of `embedding.py`) | empty file; can be deleted without effect |
| `hermite.py` | 59 | Hermite polynomial helpers — `hermite_poly`, `vibration_wave_function`, `dot` (FC overlap integral) | dormant; overlaps conceptually with `carriercapture/capture_rate.py` which computes FC overlaps via sparse-grid Schrödigner solver |
| `nu.py` | 11 | Outer product helper `direct_mul(a, b)` | dormant; `np.outer(a, b)` is the standard equivalent |
| `un.py` | 3 | Was a print-on-import snippet; now empty placeholder | dormant; no functions; safe to delete |
| `constants.py` | 2 | `h_bar` (eV·s) and `kB` (eV·K⁻¹) constants | dormant; **active code uses `from scipy import constants`** in `photoluminescence.py`, so this is unused. If revived, prefer extending `pyphotonics/constants.py` instead. |

## Planned future directions

If you (or the original author) start work on one of these, the suggested
starting point is:

- **Jahn-Teller physics** → extend `jt.py` / `jtsoc.py`. Replace the
  direct-`sympy` symbolic polynomial manipulation with numeric approaches
  if `sympy` becomes a deployment problem.
- **Defect embedding workflows** → extend `embedding.py`. Consider
  migrating to `pymatgen`'s `Structure` API rather than raw `Poscar`
  reads.
- **Generalized Huang-Rhys / harmonic analysis** → rather than reviving
  `ht.py` / `hermite.py` in isolation, integrate with the active
  `Photoluminescence` class (e.g., add a `mode_resolved_hr()` method
  that returns per-mode HR data without changing the existing API).
- **1D Schrödigner / FC overlap** → already solved in
  `carriercapture/`. Use that.

## Active grep-verified dependencies

| Module | Imports from this package | External |
|--------|---------------------------|----------|
| `jt.py` | `photonics2.ht`, `photonics2.hermite` | `numpy`, `matplotlib`, `sympy`, `joblib`, `multiprocessing`, `concurrent.futures` |
| `jtsoc.py` | `photonics2.xyz` | `numpy`, `matplotlib`, `sympy`, `joblib`, `multiprocessing`, `concurrent.futures` |
| `embedding.py` | (none) | `pymatgen`, `numpy` |
| `ht.py` | (none) | `numpy`, `matplotlib` |
| `schrodinger.py` | (none) | `scipy.linalg`, `matplotlib`, `pandas` |
| `hermite.py` | (none) | `numpy`, `matplotlib` |
| `nu.py` | (none) | `numpy` |
| `un.py` | (none) | `numpy` |
| `constants.py` | (none) | (none) |
| `ebedding.py` | (none) | (none) |

## How to verify (commands)

```bash
cd qqs/lineshape/src && \
python -c "from photonics2 import jt, jtsoc, ht, schrodinger, embedding, ebedding, hermite, nu, un, constants; print('all 10 dormant modules import cleanly')"
```

Expected: prints "all 10 dormant modules import cleanly" with no
exception, no module-level side effects.

## Recovery

If any of these modules is later deleted, recovery is via git:

```bash
git log --oneline -- qqs/lineshape/src/photonics2/<file>.py
git checkout <commit> -- qqs/lineshape/src/photonics2/<file>.py
```

Last commits touching this directory (as of 2026-07-11):
- `f265758` chore: delete lineshape_new_ref/ and simplify run_compare.py
- `69c3b30` feat: unify qqs/lineshape photoluminescence to vectorized implementation
- (pre-merge commits added the dormant modules initially)