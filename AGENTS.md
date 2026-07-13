# Agent instructions

> Entrypoint only. Prefer repo docs and code over pretraining.  
> Adapters (`CLAUDE.md`) must not carry a second full rule body.  
> **Human start:** [`README.md`](README.md).

## Precedence

1. User’s current explicit message  
2. This file  
3. Linked docs (`docs/agent-conventions.md`, package code)

## Always-on

- **What this is**: post-DFT package for defect **photoluminescence** (`pyphotonics/`) and **carrier capture** (`carriercapture/`). Upstream: Sherif Abdulkader Tawfik — GPL-3.0.  
- **Not** a VASP runner; consumes VASP/phonopy outputs.  
- **Install**: legacy `pip install -e .` (`setup.py`); Python ≥3.10. Console scripts: `pyphotonics`, `pyphotonics-incar`.  
- **No pytest suite** — “tests” are demos under `test/` / `testcode/`; validate by running demos + known-good numbers.  
- **`carriercapture`**: preserve numeric parity with Julia original where applicable.  
- **`qqs/lineshape/`**: unmaintained legacy (different style). Prefer refactored packages; merge notes in `docs/superpowers/specs/` if present.  
- Heavy dep chain via `oganesson` — expect install warnings.  
- Cite CPC 2022 paper when publishing results (see README / CLAUDE).  

## Development commands

```bash
pip install -e .

# PL CLI / API
pyphotonics -cgs CONTCAR_GS -ces CONTCAR_ES -m 189 -M phonopy -r 1000
python -c "from pyphotonics import Photoluminescence; p=Photoluminescence('CONTCAR_GS','CONTCAR_ES',189,method='phonopy'); print(p.HuangRhys)"

# FERWE/FERDO helper
pyphotonics-incar /path/to/vasp_folder

# Demo scripts (not a CI harness)
cd test/photoluminscence && python diamond.py
```

## Read on demand

| When | Read first |
|------|------------|
| Architecture, conventions, known issues | [`docs/agent-conventions.md`](docs/agent-conventions.md) |
| Human how-to / citation | [`README.md`](README.md) |
| HR divergence / tools notes | `tools/HR_DIVERGENCE.md`, `tools/RESOLUTION_NOTES.md` |
| Merge design (qqs) | `docs/superpowers/specs/` (if present) |
| Planned work | `docs/` plans under `.superpowers` / CHANGELOG |

## Keep in sync

| Topic | Files |
|-------|--------|
| Agent rules | This file canonical; `CLAUDE.md` = short + `@AGENTS.md` |
| Public API | `pyphotonics/`, `setup.py` entry points |
| Human README | usage examples ↔ real CLI flags |
