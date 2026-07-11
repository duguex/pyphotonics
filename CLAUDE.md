# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

`pyphotonics` v0.2.0 — post-processing Python package that computes defect photonic properties from VASP (DFT, constrained DFT) and phonopy outputs. Two subpackages:

- **`pyphotonics/`** — `Photoluminescence` class: Huang-Rhys factor, S(ω), PL line-shape. Driven by CLI `pyphotonics` and helper CLI `pyphotonics-incar` (writes `FERWE`/`FERDO` INCAR tags from an OUTCAR).
- **`carriercapture/`** — Carrier capture coefficients: 1D PES fitting, Schrödinger eigenstates, Franck-Condon overlaps, Marcus transfer rates, Einstein mobility. Ported from Julia (preserve numeric parity with Julia original).

Author: Sherif Abdulkader Tawfik — GPL-3.0. Cite [Sherif AbdulkaderTawfik, Salvy P.Russo, *Computer Physics Communications*, 2022, 273, 108222](https://www.sciencedirect.com/science/article/pii/S0010465521003349) when publishing results from this code.

For deep architecture, file-by-file roles, conventions, and known issues see **`AGENTS.md`** in the repo root. Do not duplicate its contents here.

## Environment

- Python ≥3.10 (uses `from __future__ import annotations`).
- Legacy `setup.py`-based install — no `pyproject.toml`, no linter/formatter/test-runner/CI.
- Runtime deps: `scipy`, `numpy`, `pandas`, `matplotlib`, `pymatgen`, `oganesson` (the last pulls a heavy chain: torch, matgl, etc. — installs with warnings).

## Install

```bash
pip install -e .
```

Registers two console scripts: `pyphotonics` (PL CLI) and `pyphotonics-incar` (FERWE/FERDO from OUTCAR).

## Run

```bash
# PL line-shape via CLI
pyphotonics -cgs CONTCAR_GS -ces CONTCAR_ES -m 189 -M phonopy -r 1000

# PL line-shape via Python API
python -c "
from pyphotonics import Photoluminescence
p = Photoluminescence('CONTCAR_GS', 'CONTCAR_ES', 189, method='phonopy')
print('S =', p.HuangRhys)
"

# FERWE/FERDO helper
pyphotonics-incar /path/to/vasp_folder

# Example scripts (NOT a test runner — just demo scripts reading disk files)
cd test/photoluminscence && python diamond.py
cd test/ferwe_ferdo && pyphotonics-incar
```

There is no `pytest`, no `unittest`, no coverage tooling. "Tests" are standalone scripts in `test/` and `testcode/`. To sanity-check after a change, run the relevant demo script and diff outputs against the user's known good values; do not look for a green test suite.

## Legacy `qqs/lineshape/` directory

Unmaintained historical code from the original author. Different conventions from the rest of the repo (no type hints, bare `except:` with `print()`+`sys.exit()`, absolute intra-package imports, `multiprocessing`/`concurrent.futures` imports). Layout:

- `src/photonics2/` — 16 Python source files (legacy PL + Jahn-Teller). `photoluminescence.py` and `plott.py` are now algorithm-equivalent to `lineshape_new_ref/` (merged 2026-07-11 per `docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`).
- `src/lineshape_new_ref/` — earlier numpy-vectorized snapshot kept for A/B regression via `run_compare.py`.
- `cases/` — 8 VASP/phonopy case studies (`1`, `123`, `beta_Ag_pair`, `Cs3Cu2Br5_STE`, `CsCuAgI3_pair`, `CuCs`, `Vbr`, `zlq`).
- `output/` — regeneratable figures and data.
- `BASELINE.md` — recorded numerical baseline (HR, Δ_R, Δ_Q) for all 8 cases; use it to verify legacy refactors.
- `run_compare.py` — side-by-side runner for the current vs vectorized implementations.

```bash
cd qqs/lineshape && python src/pl.py
python run_compare.py
```

Performance note for legacy case runs: `band.yaml` files are 27–125 MB. Pass `m=masses_kg` to skip the YAML parser bottleneck — see `run_compare.py` for an example.

Do not modernize `qqs/lineshape/` as a side effect of working in the refactored `pyphotonics/` package.

## When changing the refactored code (`pyphotonics/` and `carriercapture/`)

- Match existing conventions: full type hints, `from __future__ import annotations`, NumPy-style docstrings, `pathlib.Path` + context managers for I/O, `FileNotFoundError`/`ValueError`/`TypeError` for errors, `raise ... from exc` chaining, relative intra-package imports.
- Constants live in `pyphotonics/constants` (or `carriercapture/constants`) as `SCREAMING_SNAKE`.
- No logging framework — use `print()` or file writes.
- Single-threaded: no `multiprocessing`, `concurrent.futures`, or async additions.
- Configuration is constructor params + argparse; no env vars, no config files.
- Backward-compat aliases matter: `Photoluminescence` exposes `HuangRhyes` (typo) and `numAtoms` properties and accepts legacy constructor kwargs (`exceited_state`, `numModes`, `m`). Preserve them unless the deprecation is intentional and announced.
- `method="vasp"` is not supported and `NotImplementedError`s in that path — do not paper over it.

## Cite

If work on this repo leads to a publication, cite the Computer Physics Communications paper above (see README.md).

## Plan/spec artifacts

Active planning lives under `docs/superpowers/plans/` and `docs/superpowers/specs/` (current file naming: `YYYY-MM-DD-<topic>.md`). Treat any plan in those trees as in-flight design work — read it before starting related implementation.