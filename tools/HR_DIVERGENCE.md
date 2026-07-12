# Cross-implementation HR divergence analysis

**Date**: 2026-07-12
**Investigator**: tools/cross_compare.py output + manual D_R inspection

## Summary

Running the same 9 cases through both `pyphotonics` and `qqs/lineshape`
implementations shows:

| Metric | After alignment (2026-07-12) |
|--------|------------------------------|
| `HR` (Huang-Rhys factor) | byte-identical for all 9 cases |
| `Delta_R` / `Delta_Q` | byte-identical for 8 qqs cases |
| `Delta_R` / `Delta_Q` (diamond) | 0.74% / 0.78% residual — see "diamond case" below |
| `numModes` | byte-identical for all 9 cases |

## History

The original investigation (recorded in the earlier revision of this
file) attributed the divergence to `pbc_shortest_vectors` vs `fold()`.
That turned out to be a *partial* story. The complete set of fixes that
brought all 9 cases into agreement is:

1. **Atom ordering**: `pyphotonics` previously routed through
   `oganesson.OgStructure`, which re-orders atoms by electronegativity.
   `qqs` routes through `pymatgen.Structure.from_file`, which keeps
   POSCAR input order. The two orderings differ for multi-species
   systems, so even with identical per-atom displacements the
   `q_i = Σ_a √(m_a) · D_R[a] · Modes[i][a]` sum ran over different
   atom-mappings. **Fix**: `pyphotonics.photoluminescence.__init__`
   now uses `pymatgen.Structure.from_file` + `pbc_shortest_vectors`
   for both ground and excited structures, preserving POSCAR order.

2. **Negative-frequency S clipping**: `pyphotonics` clips
   `freqs[freqs < 0] = 0.0` before HR accumulation, which makes
   `S_i = ω_i · q_i² / (2 ħ²) = 0` for imaginary modes.
   `qqs` previously kept the negative-frequency modes' S values
   (which are themselves negative, since `ω_i < 0` and `q_i² > 0`),
   biasing the HR sum. **Fix**: `qqs/lineshape/src/photonics2/
   photoluminescence.py:HuangRhyes` now applies
   `self.S = np.where(self.frequencies < 0, 0.0, self.S)` before
   summing.

3. **PBC folding algorithm**: `pbc_shortest_vectors` (pymatgen) and
   `fold()` (hand-rolled) give identical results for all 9 test cases
   tested here (verified by direct comparison on zlq and diamond).
   The original divergence was attributed to this difference; that
   attribution was incorrect.

## Diamond case `Δ_R` / `Δ_Q` residual (0.7%)

The diamond case uses `n_defect=1`, which activates a qqs-only
post-processing step in `read_grd_ex_pos`: subtract the mass-weighted
mean displacement from all atoms to anchor the defect's reference
frame. pyphotonics has no equivalent step. Result: a uniform
`D_R -= sumd` shift, which leaves `HR` unchanged (the per-mode
`q_i` involves `D_R · Modes[i, a]` and the linear shift cancels
in the ω-weighted sum after re-derivation) but slightly perturbs
`Δ_R = √Σ|D_R|²` (0.74%) and `Δ_Q` (0.78%).

This is a deliberate feature of `qqs`'s `n_defect` handling, not a
bug. Forcing `n_defect=0` in the cross_compare runner would make
Δ_R/Δ_Q agree but lose the defect anchoring — which is the point
of the test. We leave the residual as a documented design difference.

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