# Repository Guidelines

## Project Overview

**pyphotonics** (v0.2.0) — post-processing DFT code that calculates photonic properties of defects using VASP and phonopy outputs. Two subpackages:

- **`pyphotonics`** — Photoluminescence: computes Huang-Rhys factors and PL line-shapes from ground/excited state structures + phonon modes.
- **`carriercapture`** — Carrier capture coefficients: 1D PES fitting, Schrödinger eigenstates, Franck-Condon overlaps, Marcus-theory transfer rates, Einstein mobility.

Author: Sherif Abdulkader Tawfik — GPL-3.0 — <https://github.com/sheriftawfikabbas/pyphotonics>

---

## Architecture & Data Flow

```
                         ┌─────────────────────┐
                         │  CLI (cli.py)        │
                         │  pyphotonics         │
                         │  pyphotonics-incar   │
                         └──────┬──────────────┘
                                │
             ┌──────────────────┼──────────────────┐
             ▼                  ▼                  ▼
  ┌──────────────────┐  ┌──────────────┐  ┌──────────────┐
  │ Photoluminescence│  │ CLI_INCARs   │  │ ferwe_ferdo  │
  │ (photolumines-   │  │ (ferwe/ferdo │  │ (utilities/  │
  │ cence.py)        │  │  from OUTCAR)│  │  ferwe_ferdo)│
  └────────┬─────────┘  └──────────────┘  └──────────────┘
           │
     reads VASP/phonopy
           │
           ▼
  ┌──────────────────────────────────────────────┐
  │  Output: Huang-Rhys factor, S(ω), PL lineshape│
  └──────────────────────────────────────────────┘

 carriercapture/ (standalone, included in setup.py):
   constants → brooglie → potential → capture_rate → param_scan
                                      → transfer_coord → plotter
```

**Photoluminescence pipeline**:
1. Read ground/excited state structures (pymatgen) → compute Δ vector
2. Read phonon modes and frequencies (phonopy `band.yaml`)
3. Compute Huang-Rhys factor S_k per mode, total S, Δ_R, Δ_Q, IPR
4. Compute spectral function S(ω) via Gaussian broadening
5. Compute PL line-shape via Fourier transform

**Carrier capture pipeline**:
1. Parse PES data → `Potential` dataclass (fit harmonic/Morse/polynomial/spline)
2. Filter sample points to energy island around minimum
3. Solve 1D Schrödinger equation (sparse eigen-solver) → eigenstates
4. Compute Franck-Condon overlaps → capture coefficients
5. (Optional) Marcus-type transfer rates, reorganization energy, mobility

---

## Key Directories

| Path | Purpose |
|------|---------|
| `pyphotonics/` | Main package: PL line-shape computation, CLI, version |
| `pyphotonics/utilities/` | FERWE/FERDO INCAR tag helper |
| `carriercapture/` | Carrier capture coefficient subpackage |
| `test/` | Integration/demo scripts (no test framework) |
| `testcode/` | Standalone test scripts and data |

---

## Development Commands

```bash
# Run PL calculation (CLI)
pyphotonics -cgs CONTCAR_GS -ces CONTCAR_ES -m 189 -M phonopy -r 1000

# Run PL calculation (Python API)
python -c "
from pyphotonics import Photoluminescence
p = Photoluminescence(ground_state='CONTCAR_GS', excited_state='CONTCAR_ES',
                      num_modes=189, method='phonopy')
print('S =', p.HuangRhys)
"

# Run INCAR utility
pyphotonics-incar --vasp_folder ./

# Run example scripts
cd test/photoluminscence && python diamond.py

# Install (editable)
pip install -e .
```

No test runner, linter, or formatter configured.

---

## Code Conventions & Common Patterns

### Main Package (`pyphotonics/`) — Refactored
- **Typing**: Full type hints on all functions and methods (`from __future__ import annotations`)
- **Error handling**: `FileNotFoundError` for missing input files, `ValueError`/`TypeError` for invalid params
- **Naming**: Consistent `snake_case` (methods, params, locals), `PascalCase` (classes), `SCREAMING_SNAKE` (constants)
- **Constructor pattern**: Lightweight init → internal `_compute_hr()` / `_compute_s_omega()` methods
- **Docstrings**: NumPy-style on all public methods, module-level docstring
- **Backward compat**: `HuangRhyes` property (typo alias), `numAtoms` property, constructor accepts legacy param names (`exceited_state`, `numModes`, `m`)
- **Constants**: Module-level `HBAR_JS`, `HBAR_EVS`, `AMU_TO_KG`, `THZ_TO_EV`
- **File I/O**: Uses `pathlib.Path` and context managers (`with open(...) as f`)
- **CLI**: Fixed `--contcar_excited_state` argument; all argparse attributes match constructor params
- **Dead code removed**: Empty `constants.py` deleted; unused `__init__.py` replaced with proper exports

### Carrier Capture (`carriercapture/`)
- **Typing**: `from __future__ import annotations` in all files. Mix of `typing` module and modern union syntax.
- **Dataclasses**: Core types (`Potential`, `ConfCoord`, `TransferCoord`) are `@dataclass`
- **Naming**: `snake_case` functions, `PascalCase` dataclasses, `_private` helpers, `SCREAMING_SNAKE` constants
- **Error handling**: `ValueError` for preconditions, `RuntimeError` wrapping scipy failures
- **Docstrings**: Present on all modules and public functions
- **Julia parity**: Comments note Julia original behavior

---

## Important Files

| File | Role |
|------|------|
| `pyphotonics/__init__.py` | Exports `VERSION`, `Photoluminescence` |
| `pyphotonics/__main__.py` | `python -m pyphotonics` entry (relative import) |
| `pyphotonics/cli.py` | CLI classes `CLI_Photoluminescence` and `CLI_INCARs` |
| `pyphotonics/photoluminescence.py` | `Photoluminescence` class — core PL computation |
| `pyphotonics/version.py` | Reads `VERSION` file |
| `pyphotonics/utilities/ferwe_ferdo.py` | FERWE/FERDO from OUTCAR |
| `carriercapture/__init__.py` | Re-exports all 33 public symbols |
| `carriercapture/potential.py` | `Potential` dataclass, fitting, 1D Schrödinger solver |
| `carriercapture/capture_rate.py` | `ConfCoord` — Franck-Condon overlap + capture coefficient |
| `carriercapture/transfer_coord.py` | Marcus rate, reorganization energy, mobility |
| `carriercapture/brooglie.py` | Sparse-grid Schrödinger solver backend |
| `carriercapture/param_scan.py` | Parameter scanning for capture coefficient |
| `carriercapture/plotter.py` | Matplotlib plotting for potentials and capture curves |
| `setup.py` | Installs both `pyphotonics` and `carriercapture` |

---

## Runtime/Tooling Preferences

- **Python**: >=3.10 recommended (`from __future__ import annotations`)
- **Package manager**: `pip` via `setup.py`; no `pyproject.toml` or `poetry`
- **Dependencies**: `scipy`, `numpy`, `pandas`, `matplotlib`, `pymatgen`, `oganesson`
- **No linter/formatter**: No `.flake8`, `pyproject.toml`, `isort`, `black`, or `ruff` config
- **Julia port**: `carriercapture/` is a direct port from Julia — maintain functional parity with the original

---

## Testing & QA

- **No test framework** configured (no pytest, unittest).
- **Test scripts** (`test/`, `testcode/`) are standalone `python` scripts, not runnable via test runners.
- **No CI** configuration found.
- **No coverage** tooling.
- **Manual testing**: Scripts read VASP files (CONTCAR, OUTCAR) and phonopy output from disk.

---

## Known Issues

| Issue | Location | Impact |
|-------|----------|--------|
| `vasp` method stubs removed (raise `NotImplementedError`) | `photoluminescence.py` | `method="vasp"` not supported |
| `oganesson` pulls deep dependency chain (torch, matgl, etc.) | `setup.py` | Heavy install; runs with warnings |
| No `pyproject.toml` — uses legacy `setup.py` | root | Modern packaging features unavailable |
| No test framework | project | No automated testing |
