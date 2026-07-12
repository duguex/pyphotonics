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

1. **Which atom is "atom 0" differs between libraries, not the Σ
   ordering**: `pyphotonics` previously routed through
   `oganesson.OgStructure`, which re-orders atoms by electronegativity
   (Sr before Cr for a Cr/Sr system). `qqs` routes through
   `pymatgen.Structure.from_file`, which preserves POSCAR input order
   (Cr before Sr). The `Modes` array is read from `band.yaml` in
   POSCAR input order on both sides. So `Modes[mode, atom]` is the
   same set of physical atoms in both implementations, but `D_R[i]`
   refers to a different atom: pyphotonics's `D_R[0]` is Sr's
   displacement, qqs's `D_R[0]` is Cr's. The HR formula
   `q[mode] = Σ_a √(m_a) · D_R[a] · Modes[mode, a]` is internally
   self-consistent on each side (each library matches its own `D_R`
   ordering with its own `Modes` ordering), but the two libraries
   compute a *different physical sum*: pyphotonics's `q[mode]` is
   the Sr-dominated projection, qqs's is the Cr-dominated one. The
   element-weighted sum then differs in proportion to the mode's
   participation ratio on each element. **Fix**: `pyphotonics.
   photoluminescence.__init__` now uses `pymatgen.Structure.from_file`
   + `pbc_shortest_vectors` for both ground and excited structures,
   matching qqs's POSCAR order. The element the two libraries label
   as `atom 0` is now the same.

2. **Negative-frequency S clipping**: `pyphotonics` clips
   `freqs[freqs < 0] = 0.0` before HR accumulation, which makes
   `S_i = ω_i · q_i² / (2 ħ²) = 0` for imaginary modes. `qqs`
   previously kept the negative-frequency modes' S values (which are
   themselves negative, since `ω_i < 0` and `q_i² > 0`), biasing the
   HR sum. **Fix**: `qqs/lineshape/src/photonics2/photoluminescence.
   py:HuangRhyes` now applies `self.S = np.where(self.frequencies < 0,
   0.0, self.S)` before summing.

3. **PBC folding algorithm**: `pbc_shortest_vectors` (pymatgen) and
   `fold()` (hand-rolled) give identical results for all 9 test cases
   tested here (verified by direct comparison on zlq and diamond).
   The original divergence was attributed to this difference; that
   attribution was incorrect.

### Note on `Σ_a` being commutative

The HR formula `q[mode] = Σ_a √(m_a) · D_R[a] · Modes[mode, a]` is
a sum over atoms. Reordering the array indices does not change the
result. The HR *divergence* between pyphotonics and qqs was not
caused by the two libraries disagreeing about array order — it was
caused by them disagreeing about which element to *call* `atom 0`
(and 1, 2, ..., N-1). The Σ itself is order-invariant; what changes
is the *element being summed* at each index.

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

## Verification commands

```bash
# Run the cross-implementation comparison:
python tools/cross_compare.py

# Verify atom 0 is the same element on both sides (should be the
# first atom listed in POSCAR — e.g. Cr for zlq):
python -c "
import sys, tempfile
from pathlib import Path
sys.path.insert(0, '/home/duguex/pyphotonics')
sys.path.insert(0, '/home/duguex/pyphotonics/qqs/lineshape/src')
with tempfile.TemporaryDirectory() as t:
    td = Path(t)
    (td / 'CONTCAR_GS.vasp').symlink_to('/home/duguex/pyphotonics/qqs/lineshape/cases/zlq/GRD')
    (td / 'CONTCAR_ES.vasp').symlink_to('/home/duguex/pyphotonics/qqs/lineshape/cases/zlq/EXC')
    from pymatgen.core import Structure
    gs = Structure.from_file(str(td / 'CONTCAR_GS.vasp'))
    print('atom 0:', gs.species[0], 'frac_coords:', gs.frac_coords[0])
"
```

## Related artifacts

- `tools/cross_compare.py` — runner that produced the divergence table
- `pyphotonics/photoluminescence.py` — uses
  `pymatgen.Structure.from_file` + `pymatgen.core.lattice.pbc_shortest_vectors`
  (lines 125–133).
- `qqs/lineshape/src/photonics2/photoluminescence.py` — `HuangRhyes()`
  applies `self.S = np.where(self.frequencies < 0, 0.0, self.S)` before
  summing (around line 333).