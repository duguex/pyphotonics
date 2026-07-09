# Repair Plan Implementation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all 8 remaining issues in the pyphotonics codebase across 3 delivery waves.

**Architecture:** 3 independent waves delivered sequentially. Wave 1 = quick bugfixes (no risk). Wave 2 = lazy import + CLI arg (touches import chain). Wave 3 = test infrastructure (depends on stable baseline).

**Tech Stack:** Python 3.10+, pytest (Wave 3)

## Global Constraints

- Must preserve anchor values: Delta_R=0.1486349021256089, Delta_Q=0.521858306335091, S_omega peak=88.725811 at index 59
- Python >= 3.10 compatible
- No new external dependencies beyond pytest (Wave 3)
- Each issue closed by commit message `Closes #N`

---
## Wave 1: Quick Fixes

### Task 1: Fix CLI --help crash without oganesson (#2)

**Files:**
- Modify: `pyphotonics/cli.py:10`
- Modify: `pyphotonics/cli.py:111`

- [ ] **Step 1: Move import into method body**

In `pyphotonics/cli.py`, delete line 10 (`from oganesson.io.vasp import Outcar`).
Insert it inside `CLI_INCARs.__init__` before the `Outcar()` call on line 111:

```python
class CLI_INCARs:
    def __init__(self, argv: Optional[list[str]] = None) -> None:
        # ... existing setup ...
        args = parser.parse_args(argv[1:])
        from oganesson.io.vasp import Outcar  # ← moved here
        outcar = Outcar(args.vasp_folder, "OUTCAR")
```

- [ ] **Step 2: Verify**

```bash
python -c "from pyphotonics.cli import CLI_Photoluminescence; print('OK')"
```
Expected: prints "OK", no oganesson import triggered.

- [ ] **Step 3: Commit**

```bash
git add pyphotonics/cli.py
git commit -m "fix: CLI --help crashes if oganesson not installed (#2)

Closes #2"
```

### Task 2: Fix MANIFEST.in (#3)

**Files:**
- Modify: `MANIFEST.in`

- [ ] **Step 1: Replace content**

```text
recursive-include pyphotonics *.py
recursive-include carriercapture *.py
include README.md
include LICENSE.txt
include CONTRIBUTING.rst
include CHANGELOG.md
include pyphotonics/VERSION
```

- [ ] **Step 2: Verify**

```bash
python setup.py sdist 2>&1 | grep -E "(copying|adding)"
# Should show carriercapture/ files in the archive
```

- [ ] **Step 3: Commit**

```bash
git add MANIFEST.in
git commit -m "fix: add carriercapture to MANIFEST.in for source dists (#3)

Closes #3"
```

### Task 3: Improve "vasp" method error message (#4)

**Files:**
- Modify: `pyphotonics/photoluminescence.py:139-142`

- [ ] **Step 1: Update error**

Current:
```python
raise NotImplementedError(
    f"method={method!r} is not implemented; use 'phonopy'"
)
```

Replace with:
```python
raise ValueError(
    f"unsupported method={method!r}. "
    f"Supported methods: 'phonopy', 'phonopy-siesta'. "
    f"The 'vasp' method is not implemented "
    f"(only phonopy band.yaml parsing is available)."
)
```

- [ ] **Step 2: Verify diamond.py still works**

```bash
cd test/photoluminscence && python diamond.py 2>&1 | head -3
```
Expected: Delta_R, Delta_Q, HuangRhys printed as before.

- [ ] **Step 3: Commit**

```bash
git add pyphotonics/photoluminescence.py
git commit -m "fix: clearer error for unsupported 'vasp' method (#4)

Closes #4"
```

### Task 4: Relax OCC_CUT_OFF assertion (#7)

**Files:**
- Modify: `carriercapture/capture_rate.py:12`

- [ ] **Step 1: Change constant**

```python
OCC_CUT_OFF = 1e-3  # previously 1e-5
```

- [ ] **Step 2: Verify carriercapture still works**

```bash
python -c "
from carriercapture import Potential, fit_pot, solve_pot
import numpy as np
Q = np.linspace(-2, 2, 80)
p = Potential(Q_data=Q, E_data=2*Q**2, func_type='spline', nev=10)
fit_pot(p); solve_pot(p)
print('ZPE:', p.eps[0], 'OK')
"
```
Expected: prints ZPE value.

- [ ] **Step 3: Commit**

```bash
git add carriercapture/capture_rate.py
git commit -m "fix: relax OCC_CUT_OFF from 1e-5 to 1e-3 (#7)

Closes #7"
```

### Task 5: Minor cleanup — unused imports, fig leaks, egg-info (#8)

**Files:**
- Modify: `pyphotonics/photoluminescence.py:8`
- Modify: `pyphotonics/cli.py:7`
- Modify: `carriercapture/plotter.py:14-38`

- [ ] **Step 1: Remove unused `import sys` from photoluminescence.py**

Delete line 8: `import sys`

- [ ] **Step 2: Remove unused `from pathlib import Path` from cli.py**

Delete line 7: `from pathlib import Path`

- [ ] **Step 3: Fix plotter.py fig leaks**

In `plot_pot` (line 14) and `plot_cc` (line 31), add `return ax.figure` at the end
(or `return fig, ax` when a fig is created) so callers can close the figure:

```python
def plot_pot(pot, *, plot_wf=False, ax=None, label="", scale_factor=2e-2):
    if ax is None:
        fig, ax = plt.subplots()
    # ... existing plotting ...
    return ax.figure if ax is None else ax.get_figure()

def plot_cc(cc, *, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    # ... existing plotting ...
    return ax.figure if ax is None else ax.get_figure()
```

- [ ] **Step 4: Stop tracking egg-info**

```bash
git rm --cached pyphotonics.egg-info/
```

- [ ] **Step 5: Verify no breakage**

```bash
python -c "from pyphotonics import Photoluminescence; print('OK')"
python -c "from carriercapture.plotter import plot_pot; print('OK')"
```

- [ ] **Step 6: Commit**

```bash
git add pyphotonics/photoluminescence.py pyphotonics/cli.py carriercapture/plotter.py
git commit -m "chore: cleanup unused imports, fig leaks, egg-info tracking (#8)

Closes #8"
```

---
## Wave 2: Architecture Changes

### Task 6: Lazy-load oganesson for fast import (#1)

**Files:**
- Modify: `pyphotonics/__init__.py`
- Modify: `pyphotonics/photoluminescence.py:14`

- [ ] **Step 1: Remove module-level oganesson import**

In `pyphotonics/photoluminescence.py`, delete line 14: `from oganesson import OgStructure`

- [ ] **Step 2: Add lazy import inside __init__**

Inside `Photoluminescence.__init__`, before line 120 (`og = OgStructure(...)`), add:

```python
from oganesson import OgStructure
```

- [ ] **Step 3: Make pyphotonics/__init__.py lazy**

Replace the eager import with a lazy accessor:

```python
"""PyPhotonics — post-processing DFT code for photonic properties of defects."""
from .version import VERSION


def Photoluminescence(*args, **kwargs):
    from .photoluminescence import Photoluminescence as _PL
    return _PL(*args, **kwargs)


__all__ = ["VERSION", "Photoluminescence"]
```

- [ ] **Step 4: Benchmark import time**

```bash
python -c "
import time
t0 = time.time()
import pyphotonics
print(f'import pyphotonics: {time.time()-t0:.2f}s')
t0 = time.time()
p = pyphotonics.Photoluminescence
print(f'resolve Photoluminescence: {time.time()-t0:.2f}s')
"
```
Expected: `import pyphotonics: < 0.5s` (was ~7s)

- [ ] **Step 5: Run diamond.py to verify no regression**

```bash
cd test/photoluminscence && python diamond.py 2>&1 | head -3
```
Expected: same anchor values.

- [ ] **Step 6: Commit**

```bash
git add pyphotonics/__init__.py pyphotonics/photoluminescence.py
git commit -m "perf: lazy-load oganesson for fast package import (#1)

Closes #1"
```

### Task 7: Add --path CLI argument (#5)

**Files:**
- Modify: `pyphotonics/cli.py`

- [ ] **Step 1: Add --path argument**

After the `--phonopy_path` argument (line 71), add:

```python
parser.add_argument(
    "-p",
    "--path",
    type=str,
    default="./",
    help="base directory for CONTCAR files",
)
```

- [ ] **Step 2: Pass path to Photoluminescence constructor**

On line 76, change to:

```python
pl = Photoluminescence(
    ground_state=args.contcar_ground_state,
    excited_state=args.contcar_excited_state,
    num_modes=args.num_modes,
    method=args.method,
    resolution=args.resolution,
    phonopy_path=args.phonopy_path,
    path=args.path,
)
```

- [ ] **Step 3: Verify**

```bash
python -m pyphotonics --help 2>&1 | grep -E "(-p|--path)"
```
Expected: shows `-p, --path` in help text.

- [ ] **Step 4: Commit**

```bash
git add pyphotonics/cli.py
git commit -m "feat: add --path argument for CONTCAR base directory (#5)

Closes #5"
```

---
## Wave 3: Test Infrastructure

### Task 8: Add pytest smoke tests (#6)

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/test_photoluminescence.py`
- Create: `tests/test_carriercapture.py`
- Create: `.github/workflows/ci.yml`
- Modify: `setup.py` (add pytest extra)

- [ ] **Step 1: Create pytest config**

`tests/conftest.py`:
```python
"""Shared test configuration."""
from pathlib import Path
import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def photoluminescence_data_dir() -> Path:
    return DATA_DIR
```

- [ ] **Step 2: Copy test data**

```bash
mkdir -p tests/data/phonopy
cp test/photoluminscence/CONTCAR_GS tests/data/
cp test/photoluminscence/CONTCAR_ES tests/data/
cp test/photoluminscence/phonopy/band.yaml tests/data/phonopy/
```

- [ ] **Step 3: Create PL smoke test**

`tests/test_photoluminescence.py`:
```python
"""Smoke tests for photoluminescence module."""
from pathlib import Path
import numpy as np
import pytest

from pyphotonics import Photoluminescence


def test_anchor_values():
    """Verify anchor values from diamond example."""
    p = Photoluminescence(
        ground_state="CONTCAR_GS",
        excited_state="CONTCAR_ES",
        num_modes=189,
        method="phonopy",
        phonopy_path=str(Path(__file__).parent / "data" / "phonopy"),
        path=str(Path(__file__).parent / "data"),
    )
    assert abs(p.Delta_R - 0.1486349021256089) < 1e-15
    assert abs(p.Delta_Q - 0.521858306335091) < 1e-15
    assert abs(p.HuangRhys - 2.187786) < 1e-3


def test_pl_function():
    """Verify PL function returns expected shapes."""
    p = Photoluminescence(
        ground_state="CONTCAR_GS",
        excited_state="CONTCAR_ES",
        num_modes=189,
        method="phonopy",
        phonopy_path=str(Path(__file__).parent / "data" / "phonopy"),
        path=str(Path(__file__).parent / "data"),
    )
    A, I = p.PL(2, 2, 1.95)
    assert isinstance(A, np.ndarray)
    assert isinstance(I, np.ndarray)
    assert len(A) == 5000
    assert len(I) == 5000
```

- [ ] **Step 4: Create carriercapture smoke test**

`tests/test_carriercapture.py`:
```python
"""Smoke tests for carriercapture module."""
import numpy as np
import pytest

from carriercapture import (
    Potential, fit_pot, solve_pot, find_crossing,
    ConfCoord, calc_overlap,
)


def test_potential_fit_and_solve():
    Q = np.linspace(-2, 2, 80)
    pot = Potential(Q_data=Q, E_data=2.0 * Q**2, func_type="spline", nev=10)
    fit_pot(pot)
    solve_pot(pot)
    assert len(pot.eps) == 10
    assert pot.eps[0] > 0  # positive ZPE


def test_find_crossing():
    Q = np.linspace(-2, 2, 80)
    gs = Potential(Q_data=Q, E_data=2.0 * Q**2, func_type="spline", nev=5)
    ex = Potential(Q_data=Q, E_data=0.8*(Q - 0.8)**2 + 2.0, func_type="spline", nev=5)
    fit_pot(gs); fit_pot(ex)
    solve_pot(gs); solve_pot(ex)
    q0, e0 = find_crossing(gs, ex)
    assert q0 > 0
    assert e0 > 1.0
```

- [ ] **Step 5: Add dev dependency**

In `setup.py`, add after `install_requires`:
```python
    extras_require={
        "dev": ["pytest"],
    },
```

- [ ] **Step 6: Create CI workflow**

`.github/workflows/ci.yml`:
```yaml
name: CI
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - run: pip install -e ".[dev]"
      - run: pytest -v
```

- [ ] **Step 7: Run tests**

```bash
pip install -e ".[dev]" 2>&1 | tail -3
pytest -v 2>&1
```
Expected: tests pass.

- [ ] **Step 8: Commit**

```bash
git add tests/ .github/ setup.py
git commit -m "test: add pytest smoke tests and CI workflow (#6)

Closes #6"
```

---
## Push All Waves

- [ ] **Step: Push to fork**

```bash
git push fork master
```
