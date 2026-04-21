from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

from .constants import HBAR_EV_S, KB_EV_K, TWO_PI
from .potential import Potential

OCC_CUT_OFF = 1e-5  # as in Julia


@dataclass
class ConfCoord:
    """Configuration coordinate: two potentials + coupling for capture coefficient."""

    name: str = ""
    V1: Potential = field(default_factory=Potential)
    V2: Potential = field(default_factory=Potential)
    W: float = math.inf
    g: int = 1

    eps_matrix: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    overlap_matrix: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    delta_matrix: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))

    temperature: np.ndarray = field(default_factory=lambda: np.array([]))
    capt_coeff: np.ndarray = field(default_factory=lambda: np.array([]))
    partial_capt_coeff: np.ndarray = field(default_factory=lambda: np.zeros((0, 0, 0)))


def calc_overlap(
    cc: ConfCoord,
    *,
    cut_off: float = 0.25,
    sigma: float = 0.025,
    Q0: float,
) -> None:
    """Compute vibrational overlap matrix and Gaussian delta approximation."""
    Q = cc.V1.Q
    if Q.size == 0:
        raise ValueError("cc.V1.Q is empty (fit potential first)")
    dL = (float(Q.max()) - float(Q.min())) / float(len(Q))

    n1 = cc.V1.eps.size
    n2 = cc.V2.eps.size
    cc.overlap_matrix = np.zeros((n1, n2), dtype=float)
    cc.delta_matrix = np.zeros((n1, n2), dtype=float)

    for i in range(n1):
        for j in range(n2):
            dE = abs(float(cc.V1.eps[i]) - float(cc.V2.eps[j]))
            if dE < cut_off:
                integrand = cc.V1.chi[i, :] * (cc.V1.Q - Q0) * cc.V2.chi[j, :]
                overlap = float(np.sum(integrand) * dL)
                cc.overlap_matrix[i, j] = overlap
                cc.delta_matrix[i, j] = math.exp(-(dE / sigma) ** 2 / 2.0) / (sigma * math.sqrt(2.0 * math.pi))


def calc_capt_coeff(cc: ConfCoord, V: float, temperature: Sequence[float]) -> None:
    """Compute capture coefficient vs temperature (cm^3/s).

    Parameters:
      - V: volume where coupling W is computed (cm^3)
      - temperature: iterable of temperatures in K
    """
    T = np.asarray(list(temperature), dtype=float)
    beta = 1.0 / (KB_EV_K * T)

    n1 = cc.V1.eps.size
    n2 = cc.V2.eps.size
    partial = np.zeros((n1, n2, T.size), dtype=float)

    # Partition function Z(T) over initial vibrational states
    Z = np.zeros(T.size, dtype=float)
    for eps in cc.V1.eps:
        Z += np.exp(-beta * float(eps))

    for i in range(n1):
        for j in range(n2):
            eps = float(cc.V1.eps[i])
            overlap = float(cc.overlap_matrix[i, j])
            delta = float(cc.delta_matrix[i, j])
            occ = np.exp(-beta * eps) / Z
            partial[i, j, :] = occ * (overlap ** 2) * delta

    # Prefactor before summation (Julia: V*2π/ħ*cc.g*cc.W^2)
    partial *= (V * TWO_PI / HBAR_EV_S * cc.g * (cc.W ** 2))

    # replace NaNs with 0
    partial = np.nan_to_num(partial, nan=0.0, posinf=0.0, neginf=0.0)

    # occupation of highest eigenvalue at highest temperature
    occ_high = float(np.exp(-beta[-1] * float(cc.V1.eps[-1])) / Z[-1])
    if math.isnan(occ_high):
        occ_high = 0.0
    if not (occ_high < OCC_CUT_OFF):
        raise AssertionError(f"occ(eps_max, T_max)={occ_high} should be less than {OCC_CUT_OFF}")

    capt = partial.sum(axis=(0, 1))
    # Julia replaces 0 -> 1E-127 to avoid log(0) downstream
    capt[capt == 0.0] = 1e-127

    cc.partial_capt_coeff = partial
    cc.capt_coeff = capt
    cc.temperature = T
