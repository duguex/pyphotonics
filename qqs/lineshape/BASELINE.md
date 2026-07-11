# qqs/lineshape — Numerical Baseline

Recorded 2026-07-11 at commit `d43f491`.  
Used for regression testing when refactoring either version.

## Results

All 8 cases run with `resolution=500`, `n_defect=0`, fast parser path (`m=masses_kg` passed explicitly).

| Case | Version | Huang-Rhys S | Δ_R | Δ_Q | skipmodes | num_modes |
|------|---------|-------------|-----|-----|-----------|-----------|
| Cs3Cu2Br5_STE | A (current) | 81.4196 | 1.0110 | 9.4562 | 52 | 240 |
| | B (old vec) | 81.4196 | 1.0110 | 9.4562 | 2 | 240 |
| CsCuAgI3_pair | A | 147.0123 | — | 17.1248 | 77 | — |
| | B | 147.0123 | — | 17.1248 | 15 | — |
| beta_Ag_pair | A | 18.4019 | — | 5.2283 | 48 | — |
| | B | 18.4019 | — | 5.2283 | 7 | — |
| CuCs | A | 237.8647 | — | 13.6051 | 58 | — |
| | B | 237.8647 | — | 13.6051 | 8 | — |
| Vbr | A | **-2.7232** | — | 9.5390 | 57 | — |
| | B | **-2.7232** | — | 9.5390 | 4 | — |
| zlq | A | 4.9264 | — | 1.4360 | 5 | — |
| | B | 4.9264 | — | 1.4360 | 2 | — |
| 1 | A | 8.6617 | — | 1.3707 | 2 | — |
| | B | 8.6617 | — | 1.3707 | 2 | — |
| 123 | A | 3.0521 | — | 0.6177 | 2 | — |
| | B | 3.0521 | — | 0.6177 | 2 | — |

## Notes

- **HR, Δ_R, Δ_Q identical** between versions for all 8 cases.
- **skipmodes identical** between A and B after the 2026-07-11 merge (both use threshold `freq ≤ 0.0 eV`).
- **Vbr HR is negative** (−2.72). Likely an input data issue (GS/ES ordering or convergence), not a code bug.
- **Recorded commit**: see git log; **pre-merge commit** was `d43f491`.
- Run via: `cd qqs/lineshape && python run_compare.py`
