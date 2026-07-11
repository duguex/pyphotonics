# qqs/lineshape/ Reorganization Design

**Date**: 2026-07-10
**Status**: Approved

## Goal

Reorganize the legacy `qqs/lineshape/` directory — currently a flat mix of source code, case study input data, and generated output files — into a clean three-layer structure.

## Constraints

- Preserve all files (no deletions except obvious artifacts)
- Preserve `photonics2/.git/` history
- Renaming case directories only where necessary
- Minimal changes to Python source paths

## Target Structure

```
qqs/lineshape/
├── src/                           ← all source code
│   ├── pl.py                      ← main driver script
│   ├── phonon_struct_op.py        ← phonon structure ops script
│   ├── INCAR                      ← parameter config for pl.py
│   └── photonics2/                ← subpackage (keeps its .git/)
│       ├── __init__.py
│       ├── photoluminescence.py
│       ├── plott.py
│       ├── jt.py                  ← Jahn-Teller solver
│       ├── jtsoc.py               ← JT + SOC solver
│       ├── ht.py                  ← FC factor via Hamiltonian diagonalization
│       ├── hermite.py             ← Hermite polynomial FC integrals
│       ├── configuration_coordinate.py
│       ├── embedding.py
│       ├── ebedding.py            ← (typo variant, keep)
│       ├── schrodinger.py
│       ├── xyz.py
│       ├── constants.py
│       ├── nu.py / un.py
│       └── .plott.py.swp          ← DELETE (vim swap)
├── cases/                         ← case study input data
│   ├── 1/                         ← 240 XSF modes + 2 POSCAR + band.yaml
│   ├── 123/                       ← MAO system
│   ├── Vbr/
│   ├── beta_Ag_pair/
│   ├── CuCs/
│   ├── Cs3Cu2Br5_STE/
│   ├── CsCuAgI3_pair/
│   └── zlq/
└── output/                        ← regeneratable computation outputs
    ├── data/                      ← .data files + mode_eigenvector_data
    └── figs/                      ← .png files
```

## File Inventory and Mapping

### Source moves

| From | To |
|------|----|
| `qqs/lineshape/pl.py` | `qqs/lineshape/src/pl.py` |
| `qqs/lineshape/phonon_struct_op.py` | `qqs/lineshape/src/phonon_struct_op.py` |
| `qqs/lineshape/INCAR` | `qqs/lineshape/src/INCAR` |
| `qqs/lineshape/photonics2/` | `qqs/lineshape/src/photonics2/` |

### Case directory moves (8 directories)

| From | To |
|------|----|
| `qqs/lineshape/1/` | `qqs/lineshape/cases/1/` |
| `qqs/lineshape/123/` | `qqs/lineshape/cases/123/` |
| `qqs/lineshape/Vbr/` | `qqs/lineshape/cases/Vbr/` |
| `qqs/lineshape/beta_Ag_pair/` | `qqs/lineshape/cases/beta_Ag_pair/` |
| `qqs/lineshape/CuCs/` | `qqs/lineshape/cases/CuCs/` |
| `qqs/lineshape/Cs3Cu2Br5_STE/` | `qqs/lineshape/cases/Cs3Cu2Br5_STE/` |
| `qqs/lineshape/CsCuAgI3_pair/` | `qqs/lineshape/cases/CsCuAgI3_pair/` |
| `qqs/lineshape/zlq/` | `qqs/lineshape/cases/zlq/` |

### Output file moves (.data → output/data/, .png → output/figs/)

**Data files** (13):
`AeV.data`, `C_omega.data`, `D(e-g).data`, `PLev.data`, `main_modes.data`, `mode_eigenvector_data`, `modes.data`, `partial.HuangRhyes.data`, plus `photonics2/Et.data`

**Figure files** (6):
`AeV.png`, `PLev.png`, `PLnm.png`, `Phon.png`, `Shw.png`, `Sk.png`

### Deletions

- `photonics2/.plott.py.swp` — vim swap artifact

## Path Fixes

`pl.py` contains hardcoded paths relative to CWD that point to case data:

- `./zlq/band.yaml` → `../cases/zlq/band.yaml`
- `./zlq/GS` → `../cases/zlq/GS`
- `./zlq/ES` → `../cases/zlq/ES`

The `./INCAR` default path also breaks when `pl.py` moves to `src/`. Fix by resolving INCAR relative to the script's directory:

```python
# from:  incar_path = "./INCAR"
# to:
import os
incar_path = os.path.join(os.path.dirname(__file__), "INCAR")
```

These paths are only defaults — values read from the INCAR file override them.

## Non-goals

- No renaming of GS/GRD/ES/EXC inconsistencies inside case directories
- No modification of `photonics2/` source code beyond the swap file deletion
- No changes to case data contents or structure beyond the directory move
- No CI, test, or packaging changes

## Implementation Plan

1. Create target directories (`src/`, `cases/`, `output/data/`, `output/figs/`)
2. Move source files (`pl.py`, `phonon_struct_op.py`, `INCAR`, `photonics2/`)
3. Move case directories into `cases/`
4. Move output files into `output/`
5. Delete `.plott.py.swp`
6. Update hardcoded paths in `pl.py`
7. Verify resulting structure
