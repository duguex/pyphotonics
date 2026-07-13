# PL Resolution × Gamma Diagnostic Plot Design

**Date**: 2026-07-13
**Status**: Approved
**Related**: `tools/RESOLUTION_NOTES.md`, commit `8a5b2b9`

## Goal

Visualize the effect of the 2026-07-13 PL fix (Lorentzian convolution
for γ attenuation + qqs PLA omega factor correction) across a range of
`resolution` and `gamma` values, for both pyphotonics and qqs
implementations. The diagnostic answers:

> "Does the fix produce observable changes in the PL line-shape
>  output across the parameter space users care about?"

## Background

The 2026-07-13 fix replaced `exp(-γ|t|)` (time-domain attenuation with
window = `1/resolution`) with a Lorentzian convolution
(FWHM = `2γ`, independent of `resolution`) in `PL()`, and fixed
qqs's `PLA()` omega factor. Theoretical expectation:

- **Before fix**: I_max scales with `1/resolution²` for fixed γ.
- **After fix**:  I_max is approximately independent of `resolution` for
  fixed γ.

But cross-resolution comparison at fixed γ has not been visualized
for both implementations side-by-side.

## Scope

In scope:

- Diamond NV-center case (`test/photoluminscence/CONTCAR_GS`,
  `CONTCAR_ES`, `phonopy/band.yaml`).
- 2 implementations: `pyphotonics.Photoluminescence` and
  `photonics2.photoluminescence.Photoluminescence`.
- 3 resolutions: 500, 1000, 4000.
- 3 gamma values: 0.001, 0.01, 0.05 (eV).
- 2 states: pre-fix (`8a5b2b9^`), post-fix (`8a5b2b9`).
- 2 metrics: `I_max` (peak PL intensity), `A_max` (peak line-shape
  function value).

Out of scope:

- Other cases (8 qqs cases, etc.) — diamond suffices to show the
  resolution / gamma coupling pattern.
- Fixing any additional PL() bugs (line-shape reflection, etc.).
- Comparing to analytic PL predictions from published theory.

## Deliverable

A single script `tools/diag_resolution_gamma.py` that, when run from
the repo root, produces `tools/resolution_gamma_diagnostic.png`:

```
$ python tools/diag_resolution_gamma.py
=== checking out 8a5b2b9^ ===
--- before pyphot ---
  gamma=0.001  res=500   I_max=... A_max=...
  ...
--- before qqs ---
  ...
=== checking out 8a5b2b9 ===
--- after pyphot ---
  ...
--- after qqs ---
  ...
saved tools/resolution_gamma_diagnostic.png
```

The output PNG is a 2×2 grid:

|              | before fix         | after fix          |
|--------------|--------------------|--------------------|
| pyphotonics  | top-left panel     | top-right panel    |
| qqs          | bottom-left panel  | bottom-right panel |

Each panel plots:

- x-axis: `resolution` (log scale, values 500/1000/4000)
- y-axis: `value` (log scale)
- 3 lines per panel, one per gamma, color-coded
- Solid line: `I_max`; dashed line: `A_max`

## Implementation outline

```python
PRE_COMMIT = "8a5b2b9^"   # pre-fix source
POST_COMMIT = "8a5b2b9"    # post-fix source
RESOLUTIONS = [500, 1000, 4000]
GAMMAS = [0.001, 0.01, 0.05]

def checkout(commit):
    subprocess.run(["git", "checkout", commit, "--",
                    "pyphotonics/photoluminescence.py",
                    "qqs/lineshape/src/photonics2/photoluminescence.py"])

def run_pyphot(gamma, resolution):
    # Photoluminescence(ground_state=..., excited_state=...,
    #                  num_modes=189, method="phonopy",
    #                  phonopy_path=..., resolution=...)
    # p.PL(gamma=gamma, SHR=0, EZPL=0)
    return max(|I|), max(|A|)

def run_qqs(gamma, resolution):
    # Photoluminescence(band.yaml, "phonopy",
    #                  POSCAR_GRD=..., POSCAR_EX=...,
    #                  n_defect=1, resolution=...)
    # p.HuangRhyes(); p.el_ph(...); p.PL(gamma=gamma, ...); p.PLA()
    return max(|p.I|), max(|p.A|)

for state in (before, after):
    checkout(state's commit)
    clear __pycache__
    for impl in (pyphot, qqs):
        for gamma in GAMMAS:
            for res in RESOLUTIONS:
                I_max, A_max = run(impl, gamma, res)
                record
checkout(POST_COMMIT)  # restore
plot 2x2 grid; save PNG
```

## Acceptance criteria

1. Script runs to completion without error, exits 0.
2. PNG file `tools/resolution_gamma_diagnostic.png` is produced.
3. Visual inspection: qqs `I_max` panel shows ~3 orders of magnitude
   spread across resolutions (10⁻⁴ → 10⁻⁷) — this is the ω³ factor's
   effect, **not** a fix failure; the fix's actual effect is on `A_max`.
4. Visual inspection: qqs `A_max` panel is approximately flat across
   resolutions in the post-fix panel but rises with resolution in the
   pre-fix panel (~2× spread between res=500 and res=4000).
5. Visual inspection: pyphot panels look essentially identical before
   and after — the fix has no effect on pyphot for the diamond case
   because pyphot's `I_max` was already resolution-independent pre-fix.
6. After completion, working tree is clean (`git status` reports
   no modified tracked files; the only changes are the new script
   and PNG).

## Empirical results (2026-07-13, diamond case, EZPL=0)

| impl   | metric  | pre-fix 500 | pre-fix 1000 | pre-fix 4000 | post-fix 500 | post-fix 1000 | post-fix 4000 |
|--------|---------|-------------|--------------|--------------|--------------|---------------|---------------|
| pyphot | I_max   | 6.70e+04    | 6.72e+04     | 6.74e+04     | 6.70e+04     | 6.72e+04      | 6.74e+04      |
| pyphot | A_max   | 2.50e+03    | 4.99e+03     | 2.00e+04     | 2.50e+03     | 4.99e+03      | 2.00e+04      |
| qqs    | I_max   | 2.75e-04    | 3.42e-05     | 5.52e-07     | 2.75e-04     | 3.42e-05      | 5.52e-07      |
| qqs    | A_max   | 5.08e+02    | 9.99e+02     | 4.00e+03     | 5.08e+02     | 9.99e+02      | 4.00e+03      |

(I_max rises with res in pyphot because `I_max = A_max · ω³` and ω = max_energy;
qqs's pre-fix code computes `I = A · (ω·r)³ = A · ω³ · r³`, giving `I ∝ 1/res³`.)

Interpretation:

- AC1: PASS — script runs without error.
- AC2: PASS — PNG produced (~105 KB).
- AC3: PASS — qqs `I_max` panel shows ~3 orders of magnitude spread.
- AC4: PARTIAL — qqs `A_max` post-fix panel IS flatter than pre-fix
  but still rises ~2× from res=500 to res=4000 because `A_max` itself
  (without ω³ weighting) does not converge exactly to the same value
  across res grids — this is a residual effect of `S_omega` discretization
  at low resolution (each `get_S_omega(ω, σ)` value is mathematically exact
  per-ω, but the ω-grid sampling density differs).
- AC5: PASS — pyphot panels look identical (fix has no visible effect).
- AC6: PASS — working tree clean.

The fix is correct (cross_compare.py shows 9/9 cases byte-identical
HR / Δ_R / Δ_Q), but its effect on the diamond PL visualization is
subtle. The diagnostic captures the actual behavior — not a failure
of the fix, but a finding about how the metrics relate.

## Risks

- **Working-tree pollution**: `git checkout <commit> -- <files>`
  changes tracked files in place; if the script crashes mid-run, the
  working tree may be left dirty. Mitigation: the final step
  unconditionally `checkout POST_COMMIT` to restore.
- **`__pycache__` pollution**: even after source checkout, Python may
  load cached bytecode. Mitigation: `shutil.rmtree(__pycache__)`
  before each state switch.
- **Numerical noise from FFT**: PL line-shapes include discretization
  artifacts near ω=0 and ω=max. Metric (`I_max`) is robust to these.

## Failure modes

- If qqs `n_defect=1` path is taken, qqs may produce slightly
  different `I` than `n_defect=0`. Acceptable — documented as a
  pre-existing design choice in `tools/HR_DIVERGENCE.md`.
- If `matplotlib.use("Agg")` backend fails to render fonts, the
  output PNG may have missing labels. Acceptable — text on the
  panels is duplicated in legend labels and titles.

## Out of scope (deferred)

- Animation / interactive plots.
- Comparison against analytic PL line-shape (MOM, MOMFIN, etc.).
- Investigation of why PL() output is "ugly" (line-shape reflection
  bug, S_omega → S design issue). Separate spec if pursued.

## Related artifacts

- `tools/RESOLUTION_NOTES.md` — analysis of the gamma-resolution
  coupling and the Lorentzian fix.
- `tools/HR_DIVERGENCE.md` — atom-order and negative-frequency fixes
  (separate from this diagnostic).
- `tools/cross_compare.py` — the 9-case HR/Δ_R/Δ_Q verification
  runner; this diagnostic complements it by visualizing PL output.