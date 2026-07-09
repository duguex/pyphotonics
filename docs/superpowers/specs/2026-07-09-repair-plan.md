# Repair Plan: pyphotonics Remaining Issues

Date: 2026-07-09
Author: duguex
Status: Draft

## Overview

Fix 8 open GitHub issues on `github.com/duguex/pyphotonics`, organized into 3 delivery
waves by dependency and impact. After this plan is executed, all known quality issues
from the initial refactoring (commit `1801681`) will be resolved.

---

## Wave 1 — Quick Fixes (5 issues, ~2 min each)

### 1.1 Issue #2: CLI --help crashes without oganesson

**File:** `pyphotonics/cli.py`

**Problem:** `from oganesson.io.vasp import Outcar` at module level means
`pyphotonics --help` fails if oganesson isn't installed.

**Fix:** Move the import inside `CLI_INCARs.__init__`, right before the `Outcar()`
call.

**Verification:** `python -c "from pyphotonics.cli import CLI_Photoluminescence"` should
work without oganesson.

### 1.2 Issue #3: MANIFEST.in missing carriercapture

**File:** `MANIFEST.in`

**Problem:** Only `pyphotonics/VERSION` is included. `carriercapture/` has no entry,
so `python setup.py sdist` omits it.

**Fix:** Replace existing content with:
```
recursive-include pyphotonics *.py
recursive-include carriercapture *.py
include pyphotonics/VERSION
include README.md
include LICENSE.txt
include CONTRIBUTING.rst
include CHANGELOG.md
```

### 1.3 Issue #4: "vasp" method error message

**File:** `pyphotonics/photoluminescence.py:139-142`

**Problem:** `NotImplementedError` with a terse message. The "vasp" method is documented
everywhere but users get an unhelpful error.

**Fix:** Change to `ValueError` with explicit list of supported methods and explanation
that only phonopy band.yaml parsing is implemented.

### 1.4 Issue #7: carriercapture strict assertion

**File:** `carriercapture/capture_rate.py:12`

**Problem:** `OCC_CUT_OFF = 1e-5` is too strict for prototyping with few eigenstates.

**Fix:** Relax to `OCC_CUT_OFF = 1e-3`. This still guarantees convergence for real use
but doesn't block exploration.

### 1.5 Issue #8: Minor cleanup

- `pyphotonics/photoluminescence.py:8` — remove unused `import sys`
- `pyphotonics/cli.py:7` — remove unused `from pathlib import Path`
- `carriercapture/plotter.py` — `plot_pot` and `plot_cc` return `ax.figure`
  so callers can close figures and prevent memory leaks
- `git rm --cached pyphotonics.egg-info/` followed by `.gitignore` addition
  (already has `*.egg-info/`) to stop tracking build artifacts

---

## Wave 2 — Architecture Changes (2 issues)

### 2.1 Issue #1: Slow import (7s due to oganesson/torch)

**Problem:** Two levels of eager loading:

1. `pyphotonics/__init__.py` does `from .photoluminescence import Photoluminescence`
   at module level

2. `photoluminescence.py` does `from oganesson import OgStructure` at module level,
   triggering torch (~915MB), matgl, torch-geometric, lightning, etc.

**Fix (level 1):** `pyphotonics/__init__.py` — replace direct import with a lazy
accessor function:

```python
def Photoluminescence(*args, **kwargs):
    from .photoluminescence import Photoluminescence as _PL
    return _PL(*args, **kwargs)
```

This way `from pyphotonics import VERSION` stays instant but `from pyphotonics import
Photoluminescence` triggers the (still-slow) import only when actually used.

**Fix (level 2):** `photoluminescence.py` — move `from oganesson import OgStructure`
inside `__init__`, right before it's used:

```python
def __init__(self, ...):
    from oganesson import OgStructure
    og = OgStructure(file_name=str(gs_file))
```

This is safe because `OgStructure` is only needed at construction time.

**Performance targets:**
- `import pyphotonics`: < 0.5s (was ~7s)
- `Photoluminescence(...)` construction: ~7s unavoidable (oganesson + structure loading)

### 2.2 Issue #5: CLI missing --path argument

**File:** `pyphotonics/cli.py`

**Problem:** `Photoluminescence.__init__` accepts `path` (CONTCAR base directory) but
the CLI has no corresponding argument — users must `cd` to the data directory.

**Fix:** Add `-p`/`--path` argument with `default="./"` to `CLI_Photoluminescence`.
Pass to `Photoluminescence(path=args.path, ...)`.

---

## Wave 3 — Test Infrastructure (1 issue)

### 3.1 Issue #6: No test framework

**New files:**
- `tests/conftest.py` — shared fixtures (path to test data)
- `tests/test_photoluminescence.py` — smoke test that runs diamond example
  and asserts anchor values:
  - `Delta_R == 0.1486349021256089`
  - `Delta_Q == 0.521858306335091`
  - `HuangRhys ≈ 2.1878`
- `tests/test_carriercapture.py` — tests for Potential, solve_pot, capture_rate,
  transfer_coord public functions
- `.github/workflows/ci.yml` — GitHub Actions: pip install → pytest

**Test data:** Copy `test/photoluminscence/{CONTCAR_GS,CONTCAR_ES,phonopy/}` into
`tests/data/`.

**Dependencies:** Add `pytest` to setup.py `extras_require = {"dev": ["pytest"]}`.

---

## Implementation Order

```
Wave 1: #2 → #3 → #4 → #7 → #8  (parallel-safe, no dependencies)
          │
          ▼
Wave 2: #1 → #5                  (#5 depends on #1 being stable)
          │
          ▼
Wave 3: #6                        (depends on Wave 2 baseline)
```

Each wave is a separate commit with message format:
```
fix: <issue title> (#N)

Closes #N
```

---

## Verification

After each wave, run:

```bash
# Wave 1
python -c "from pyphotonics.cli import CLI_Photoluminescence; print('CLI import OK')"
python setup.py sdist  # verify carriercapture/ included
cd test/photoluminscence && python diamond.py  # anchor values unchanged

# Wave 2
python -c "import pyphotonics; print('fast import')"
python -m pyphotonics --help  # shows --path arg

# Wave 3
pip install -e ".[dev]" && pytest -v
```
