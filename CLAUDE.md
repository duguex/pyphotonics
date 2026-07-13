# CLAUDE.md

Claude Code loads this file every session. **Canonical rules:** [`AGENTS.md`](AGENTS.md).  
Local block = fallback if import is declined + high-frequency notes.

> **Humans:** [`README.md`](README.md).

## Local hard constraints (fallback)

- Two maintained packages: `pyphotonics/` (PL / Huang–Rhys) and `carriercapture/` (1D PES → FC → rates).  
- **No pytest** — demos in `test/` / `testcode/` only.  
- Do not “modernize” `qqs/lineshape/` unless the user asked; treat as legacy.  
- Preserve `carriercapture` numeric behavior vs Julia port.  
- Depth: `docs/agent-conventions.md`.

## Commands (quick)

```bash
pip install -e .
pyphotonics -cgs CONTCAR_GS -ces CONTCAR_ES -m 189 -M phonopy -r 1000
pyphotonics-incar /path/to/vasp_folder
```

## Cite

Sherif Abdulkader Tawfik, Salvy P. Russo, *Computer Physics Communications*, 2022, 273, 108222.

## Claude Code notes

- If prompted to approve imports, **allow `@AGENTS.md`**.

@AGENTS.md
