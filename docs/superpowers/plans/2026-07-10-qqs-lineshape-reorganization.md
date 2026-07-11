# qqs/lineshape/ Reorganization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure `qqs/lineshape/` from a flat directory into three layers: `src/` (code), `cases/` (input data), `output/` (regeneratable outputs).

**Architecture:** Pure file moves + path fixes in `pl.py`. No functional changes.

**Tech Stack:** Shell (mv), Python (path edit in pl.py)

## Global Constraints

- Preserve `photonics2/.git/` entirely
- No case data content changes — only directory moves
- Delete only `.plott.py.swp`

---

### Task 1: Create target directory structure

**Files:**
- Create: `qqs/lineshape/src/`
- Create: `qqs/lineshape/src/photonics2/`
- Create: `qqs/lineshape/cases/`
- Create: `qqs/lineshape/output/data/`
- Create: `qqs/lineshape/output/figs/`

**Interfaces:**
- Consumes: nothing
- Produces: empty target directories for subsequent move tasks

- [ ] **Step 1: Create directories**

```bash
mkdir -p qqs/lineshape/src qqs/lineshape/cases \
         qqs/lineshape/output/data qqs/lineshape/output/figs
```

---

### Task 2: Move source code to src/

**Files:**
- Move: `qqs/lineshape/pl.py` → `qqs/lineshape/src/pl.py`
- Move: `qqs/lineshape/phonon_struct_op.py` → `qqs/lineshape/src/phonon_struct_op.py`
- Move: `qqs/lineshape/INCAR` → `qqs/lineshape/src/INCAR`
- Move: `qqs/lineshape/photonics2/` → `qqs/lineshape/src/photonics2/`
- Delete: `qqs/lineshape/src/photonics2/.plott.py.swp`

**Interfaces:**
- Consumes: empty `src/` from Task 1
- Produces: all source code under `src/`

- [ ] **Step 1: Move files**

```bash
cd qqs/lineshape
mv pl.py phonon_struct_op.py INCAR src/
mv photonics2/ src/
```

- [ ] **Step 2: Delete vim swap**

```bash
rm -f src/photonics2/.plott.py.swp
```

---

### Task 3: Move case data to cases/

**Files:**
- Move: `qqs/lineshape/{1,123,Vbr,beta_Ag_pair,CuCs,Cs3Cu2Br5_STE,CsCuAgI3_pair,zlq}/` → `qqs/lineshape/cases/{...}/`

**Interfaces:**
- Consumes: empty `cases/` from Task 1
- Produces: all 8 case directories under `cases/`

- [ ] **Step 1: Move case directories**

```bash
cd qqs/lineshape
mv 1/ 123/ Vbr/ beta_Ag_pair/ CuCs/ Cs3Cu2Br5_STE/ CsCuAgI3_pair/ zlq/ cases/
```

---

### Task 4: Move output files to output/

**Files:**
- Move root `.data` + `mode_eigenvector_data` → `output/data/`
- Move root `.png` → `output/figs/`
- Move `src/photonics2/Et.data` → `output/data/`

**Interfaces:**
- Consumes: empty `output/{data,figs}/` from Task 1
- Produces: output files categorized

- [ ] **Step 1: Move data files**

```bash
cd qqs/lineshape
mv *.data mode_eigenvector_data output/data/
mv *.png output/figs/
mv src/photonics2/Et.data output/data/
```

---

### Task 5: Fix hardcoded paths in pl.py

**Files:**
- Modify: `qqs/lineshape/src/pl.py`

**Interfaces:**
- Consumes: Task 2-4 directory layout
- Produces: `pl.py` works out of the box after restructuring

- [ ] **Step 1: Fix `./zlq/` defaults**

Change line 130-132 from:
```python
    incar_params.get("BAND_YAML", "./zlq/band.yaml"),
    incar_params.get("GRD_POSCAR", "./zlq/GS"),
    incar_params.get("EXC_POSCAR", "./zlq/ES"),
```
to:
```python
    incar_params.get("BAND_YAML", "../cases/zlq/band.yaml"),
    incar_params.get("GRD_POSCAR", "../cases/zlq/GS"),
    incar_params.get("EXC_POSCAR", "../cases/zlq/ES"),
```

- [ ] **Step 2: Fix `./INCAR` default**

Change line 123 from:
```python
incar_path = "./INCAR"
```
to:
```python
import os
incar_path = os.path.join(os.path.dirname(__file__), "INCAR")
```

---

### Task 6: Verify resulting structure

**Files:**
- Verify: complete tree under `qqs/lineshape/`

- [ ] **Step 1: Verify layout matches spec**

```bash
cd qqs/lineshape
echo "=== src/ ===" && ls src/ && ls src/photonics2/*.py | head -5 && \
echo "=== cases/ ===" && ls cases/ && \
echo "=== output/data/ ===" && ls output/data/ | head -5 && \
echo "=== output/figs/ ===" && ls output/figs/ | head -5 && \
echo "=== total py: $(find . -name '*.py' | wc -l) ===" && \
echo "=== total data: $(find output -name '*.data' | wc -l) ===" && \
echo "=== total png: $(find output -name '*.png' | wc -l) ===" && \
echo "=== total cases: $(ls -d cases/*/ | wc -l) ==="
```

- [ ] **Step 2: Verify no stray files left in root**

```bash
ls qqs/lineshape/*.py qqs/lineshape/*.data qqs/lineshape/*.png 2>&1
# Expected: "No such file or directory" for each glob
```
