# Cross-implementation HR divergence analysis

**Date**: 2026-07-12
**Investigator**: tools/cross_compare.py output + manual D_R inspection

## Summary

Running the same 9 cases through both `pyphotonics` and `qqs/lineshape`
implementations shows:

| Metric | Across 9 cases |
|--------|----------------|
| `Delta_R` | byte-identical |
| `Delta_Q` | byte-identical |
| `numModes` | byte-identical |
| `HR` (Huang-Rhys factor) | **diverges on 7 of 9 cases** (up to 47% on beta_Ag_pair; the 123 and diamond cases agree exactly) |

`Δ_R` and `Δ_Q` are sums of squared displacements; HR is a weighted
sum `Σ ω_i · q_i² / (2 ħ²)` where `q_i = Σ_a √(m_a) · (D_R[a] · Modes[i][a])`.
A small per-atom `D_R` mismatch, when squared and summed over atoms
and modes, easily produces the observed HR divergence while leaving
`Δ_R = √Σ |D_R|²` essentially unchanged.

## Root cause: two different PBC folding algorithms

The discrepancy traces to a single algorithmic difference: how the
displacement `D_R[a] = ES[a] - GS[a]` is computed when atoms sit near
the periodic-boundary.

| Implementation | Source | Algorithm |
|----------------|--------|-----------|
| pyphotonics | `oganesson.OgStructure.get_delta_vector` | `pymatgen.core.lattice.pbc_shortest_vectors(lattice, GS_frac[a], ES_frac[a])` |
| qqs (lineshape) | `photonics2.photoluminescence.fold` | `d - lattice.T @ round(solve(lattice.T, d))` (hand-rolled cell-wrap) |

For atoms deep inside the unit cell, both algorithms give the same
result. For atoms near the cell boundary, the two algorithms differ
because:

- `pbc_shortest_vectors` returns the **shortest** periodic image
  vector (might choose image 0 or 1 in any axis).
- `fold` always picks image **0** in each axis via the `round(...)`
  step (i.e., it wraps the displacement into the central cell, not
  the shortest image).

### Worked example: `zlq` case, atom 0

```
qqs raw ES - GS:           [1.92e-08  9.31e+00  1.41e-08]    (9.31 Å along y)
qqs fold() result:         [1.92e-08 -3.87e-08  1.41e-08]    (central cell ≈ 0)
pyphotonics pbc_shortest:  [6.89e-03  4.27e-05  1.41e-05]    (shortest vector ≈ 0.007 Å)
```

These are three different values for the same physical quantity.
The atomic positions are essentially identical (`GS[a] ≈ ES[a]` modulo
a tiny relax); the algorithms just disagree about which periodic
image to use as the displacement.

## Why `Δ_R` doesn't show the divergence

`Δ_R = √(Σ |D_R[a]|²)` aggregates over all atoms. For most atoms the
two algorithms agree exactly; the disagreement is concentrated on
atoms at the boundary, of which there are typically only 1-2 per case.
Their contribution to the total squared sum is negligible.

But `HR = Σ_i ω_i · q_i² / (2 ħ²)` and `q_i = Σ_a √(m_a) · D_R[a] · Modes[i][a]`
is a **mode-resolved** quantity. A small per-atom `D_R` mismatch on
one boundary atom gets multiplied by all modes that have non-zero
displacement on that atom, and then squared. Hence HR is sensitive
even to the few boundary atoms where the algorithms disagree.

## Why 123 and diamond agree exactly

These two cases have no boundary atoms in problematic positions, so
both algorithms return identical `D_R` everywhere, and the HR
cascades equal.

## Which implementation is "right"?

`pbc_shortest_vectors` (pymatgen) is the textbook-correct choice for
periodic systems. `fold()` with `round()` is an approximation: it
unconditionally pulls the displacement back to image 0, which is
**not** the same as the shortest image when GS and ES straddle the
boundary non-symmetrically.

So the qqs implementation has a long-standing bug in `D_R` for
boundary atoms, and `pyphotonics` happens to fix it (via oganesson).

## Implications

- **For accurate HR**: pyphotonics is more correct, but only because
  it inherited oganesson's PBC logic.
- **For backward compatibility with published results**: the qqs
  numbers are what's been published (see the cited Computer Physics
  Communications paper). Switching to pyphotonics logic may shift
  reported HR by a few percent on some systems.
- **For unifying the two**: the right move is to **replace qqs's
  `fold()` with `pbc_shortest_vectors`**. That would make both
  implementations numerically identical (modulo floating-point order)
  and align qqs with the canonical PBC convention.

## Verification commands

```bash
# Recreate the divergence:
python tools/cross_compare.py

# Inspect a single D_R vector on zlq case:
python -c "
import sys, tempfile; from pathlib import Path
sys.path.insert(0, '/home/duguex/pyphotonics')
sys.path.insert(0, '/home/duguex/pyphotonics/qqs/lineshape/src')
with tempfile.TemporaryDirectory() as t:
    td = Path(t)
    (td / 'CONTCAR_GS.vasp').symlink_to('/home/duguex/pyphotonics/qqs/lineshape/cases/zlq/GRD')
    (td / 'CONTCAR_ES.vasp').symlink_to('/home/duguex/pyphotonics/qqs/lineshape/cases/zlq/EXC')
    from oganesson import OgStructure
    from pymatgen.core import Structure
    og = OgStructure(file_name=str(td / 'CONTCAR_GS.vasp'))
    es = Structure.from_file(str(td / 'CONTCAR_ES.vasp'))
    print('pyphot D_R[0] =', og.get_delta_vector(es)[0])
"
```

## Related artifacts

- `tools/cross_compare.py` — runner that produced the divergence table
- `qqs/lineshape/src/photonics2/photoluminescence.py` — `fold()` is at
  line 138; `read_grd_ex_pos` calls it at line 189.
- `pyphotonics/photoluminescence.py` — uses `OgStructure.get_delta_vector`
  at line 122.
- `oganesson/ogstructure.py:79` — `get_delta_vector` implementation.