# Resolution and gamma: numerical coupling in PL line-shape

**Date**: 2026-07-12
**Investigator**: tools/cross_compare.py + manual tests on zlq/diamond

## TL;DR

- **HR (Huang-Rhys factor)** is independent of `resolution` in both
  pyphotonics and qqs — verified on diamond and zlq cases.
- **S(ω) values at any single energy ω are independent of `resolution`**
  in pyphotonics (each is computed analytically by summing over modes).
- **PL line-shape (I) is sensitive to `resolution`** in **both**
  libraries, because the FFT / inverse-FFT and the gamma
  attenuation window are coupled to the energy grid.

## Verified: HR is resolution-independent

Tested with `tools/cross_compare.py`-style invocations at
`resolution ∈ {500, 1000, 4000}`:

| case | res=500 | res=1000 | res=4000 |
|------|---------|----------|----------|
| diamond | 2.187786 | 2.187786 | 2.187786 |
| zlq | 4.926363 | 4.926363 | 4.926363 |

`Δ_R` and `Δ_Q` likewise byte-identical. The reason: `_compute_hr`
sums over modes directly, never touches an energy grid.

## Verified: S(ω) point values are resolution-independent

`pyphotonics.get_S_omega(omega, sigma)` sums a Gaussian over each
mode's contribution:

```python
total = sum(s_k * exp(-(omega - freq_k)**2 / (2 sigma**2)) for each mode)
```

The sum is over modes (3N = 468 in zlq), not over an energy grid.
Tested:

```
S_omega[0]  res=500: 10.09629    res=1000: 10.09629    res=4000: 10.09629
S_omega[-1] res=500: 0          res=1000: 0          res=4000: 0
```

Identical to all printed digits. The line-shape *function* is
unchanged; only the *grid on which it's sampled* differs.

## Issue: PL line-shape couples to `resolution` via `gamma`

The PL() method applies a Gaussian attenuation
`exp(-γ |t|)` in the time domain before the inverse FFT. The `t`
axis is:

```python
# qqs/lineshape/src/photonics2/photoluminescence.py PL()
r = 1 / resolution
t = r * (np.arange(n) - n / 2)         # t spans [-1/(2r), +1/(2r)]
Gt = exp(St + Ct + ... - SHR - gamma * |t|)
A = fft.ifft(Gt) * n / (max_energy - min_energy)
```

The attenuation window is exactly `[-1/(2·resolution), +1/(2·resolution)]`.
**Higher `resolution` → narrower window → sharper lineshape**. This is a
*real coupling*, not a bug, but users need to know:

- I_max for zlq at γ=0.01 eV: 5.04e-5 (res=500) → 7.90e-6 (res=1000) → 4.94e-7 (res=4000). Scales as 1/resolution².
- The integral (area under I) is approximately preserved when
  resolution doubles (FFT preserves total area), but the peak height
  shrinks quadratically with resolution.

**pyphotonics has the same coupling** in its `PL()` method:

```python
# pyphotonics/photoluminescence.py PL()
Gt[i] = G[i] * np.exp(-gamma * |t|)
t = r * (i - n / 2)
```

`gamma * |t|` is in eV·steps. Same 1/resolution window.

## Practical guidance

For PL line-shape comparison between two implementations:

- Keep `resolution` fixed across both runs (cross_compare.py does this).
- Keep `gamma` fixed.
- The integration window `[-1/(2r), +1/(2r)]` must comfortably exceed
  `1/γ`. If `γ` is very small (long-lived ZPL) and `resolution` is
  high (narrow window), `Gt` is barely attenuated and the line-shape
  looks unphysically sharp.

For Huang-Rhys comparisons, `resolution` doesn't matter.

## Related code

- `qqs/lineshape/src/photonics2/photoluminescence.py:232-263` — PL()
  method with γ attenuation.
- `pyphotonics/photoluminescence.py:335-381` — PL() method with same
  coupling.
- `tools/cross_compare.py` — runner; uses `resolution=500` for the 8
  qqs cases and `resolution=1000` for the diamond case.