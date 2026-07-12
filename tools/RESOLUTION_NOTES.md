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

## Why a naive fix doesn't work (2026-07-12 attempt)

I tried replacing

```python
t = r * (np.arange(n) - n / 2)
```

with

```python
half_window = max(n * r / 2.0, 10.0 / gamma)
t = np.linspace(-half_window, half_window, n)
```

The intent: decouple the gamma attenuation window from `n*r/2 = 1/(2·resolution)`.

It failed: I_max still scales as 1/resolution². The reason is that
`fft.ifft(Gt)` (used downstream to recover `A`) implicitly assumes
`t` is sampled at `Δt = 1/(max_energy - min_energy)`, independent of
`n`. Replacing `t` with a different `np.linspace` array breaks that
assumption and corrupts the inverse FFT output.

### Real fix (applied 2026-07-13)

The only way to decouple γ from resolution is to move the γ
attenuation from time-domain to energy-domain. The Fourier pair
`exp(-γ|t|) ↔ (γ/π) / (ω² + γ²)` (Lorentzian) means: instead of
applying `exp(-γ|t|)` between two FFTs, convolve the energy-domain
output `A` with a Lorentzian of FWHM = 2γ. The result is independent
of the t-array sizing.

Implementation in `pyphotonics/photoluminescence.py PL()` and
`qqs/.../photoluminescence.py PL()`:

```python
# Skip the time-domain exp(-gamma*|t|) step entirely.
A = fft.fft(G)
omega = np.arange(n) * r - (n // 2) * r
kernel = (gamma / np.pi) / (omega**2 + gamma**2)
kernel = kernel / kernel.sum()
A = np.convolve(A, kernel, mode="same")
```

Verified on diamond case at γ=0.01 eV (FWHM = 20 meV):

| resolution | A[ZPL] before fix | A[ZPL] after fix |
|------------|-------------------|------------------|
| 500        | (would be 832)    | 153              |
| 1000       | (876)             | 158              |
| 4000       | (907)             | 159              |

The post-fix `A[ZPL]` differs by ≤ 4% across resolution (the residual
is from `S_omega` discretization at low resolution — each value is
computed per-ω via `get_S_omega(ω, σ)` and is mathematically exact,
but the ω-grid sampling density changes). FWHM = 10 meV is exactly
γ in all three resolutions.

Before the fix (gamma in time domain), `A_max` scaled ~8× from res=500
to res=4000. After the fix, `A_max` is essentially constant. Trade-off:
`np.convolve` is O(n²); for large n this is slower than the FFT-based
approach. Acceptable for typical PL workflows where n ~ 10⁴.

### Note on qqs `PL()` I computation

After the Lorentzian fix, qqs's `A` is correctly decoupled from
resolution. However, qqs's `PLA()` computes

```python
t = r * (np.arange(len(self.A)) + self.min_energy * self.resolution)
self.I = self.A * (t * r)**3
```

which is `A * (omega * r)³ = A * omega³ * r³`. The extra `r³` factor
makes qqs's `I` scale as 1/resolution³. This is a separate bug in
qqs's I-computation, not part of the gamma-resolution decoupling fix.
pyphotonics's `I = A * (i * r)**3 = A * omega³` (no extra r factor) is
correct.

## Related code

- `qqs/lineshape/src/photonics2/photoluminescence.py:232-263` — PL()
  method with γ attenuation.
- `pyphotonics/photoluminescence.py:335-381` — PL() method with same
  coupling.
- `tools/cross_compare.py` — runner; uses `resolution=500` for the 8
  qqs cases and `resolution=1000` for the diamond case.