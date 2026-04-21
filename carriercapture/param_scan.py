"""Parameter scan helpers (ported from ParamScan.jl)."""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

from .constants import AMU_EV_C2, HBARC_EV_M
from .potential import Potential, harmonic, morse, solve_pot
from .capture_rate import ConfCoord, calc_overlap, calc_capt_coeff


def get_harmonic_q_m(hw: float, E0: float) -> float:
    """ΔQ for activationless Marcus regime for two harmonic PES."""
    a = AMU_EV_C2 / 2.0 * (hw / HBARC_EV_M / 1e10) ** 2
    return math.sqrt(E0 / a)


def get_morse_q_m(a: float, b: float, E0: float) -> float:
    """ΔQ for activationless Marcus regime for two Morse PES."""
    # Julia: Q_m = (1/b)*log(1-(E0/a)^0.5)
    return (1.0 / b) * math.log(1.0 - math.sqrt(E0 / a))


def fit_harmonic_params(hw_i: float, hw_f: float, dQ: float, dE: float) -> Tuple[Potential, Potential]:
    """Construct and solve two harmonic potentials (excited/relaxed)."""
    nev_excited = 180
    nf_min = (1.0 / hw_f) * (((nev_excited + 0.5) * hw_i) + 2.0 * dE) - 0.5
    nev_relaxed = int(math.floor(nf_min + 1.0))

    Q = np.linspace(-20.0 - dQ, 20.0 + dQ, 5000)

    potf = Potential(name="Relaxed state")
    potf.Q0 = dQ
    potf.E0 = 0.0
    potf.nev = nev_relaxed
    potf.func = lambda x: harmonic(x, hw_f, E0=potf.E0, Q0=potf.Q0)
    potf.Q = Q
    potf.E = potf.func(Q)
    solve_pot(potf)

    poti = Potential(name="Excited state")
    poti.Q0 = 0.0
    poti.E0 = dE
    poti.nev = nev_excited
    poti.func = lambda x: harmonic(x, hw_i, E0=poti.E0, Q0=poti.Q0)
    poti.Q = Q
    poti.E = poti.func(Q)
    solve_pot(poti)

    return poti, potf


def fit_morse_params(a_i: float, a_f: float, b_i: float, b_f: float, Q0: float, E0: float) -> Tuple[Potential, Potential]:
    """Construct and solve two Morse potentials (excited/relaxed)."""
    nev_excited = 180
    nev_relaxed = 800

    Q = np.linspace(-20.0 - Q0, 20.0 + Q0, 5000)

    poti = Potential(name="Excited state")
    poti.Q0 = 0.0
    poti.E0 = E0
    poti.nev = nev_excited
    poti.func = lambda x: morse(x, np.array([a_i, b_i]), E0=poti.E0, Q0=poti.Q0)
    poti.Q = Q
    poti.E = poti.func(Q)
    solve_pot(poti)

    potf = Potential(name="Relaxed state")
    potf.Q0 = Q0
    potf.E0 = 0.0
    potf.nev = nev_relaxed
    potf.func = lambda x: morse(x, np.array([a_f, b_f]), E0=potf.E0, Q0=potf.Q0)
    potf.Q = Q
    potf.E = potf.func(Q)
    solve_pot(potf)

    return poti, potf


def get_harmonic_capture(hw_i: float, hw_f: float, dQ: float, dE: float, W: float, V: float, Tmin: float = 100.0, Tmax: float = 800.0, NT: int = 100):
    poti, potf = fit_harmonic_params(hw_i, hw_f, dQ, dE)
    cc = ConfCoord(V1=poti, V2=potf, W=W, g=1)
    # Choose Q0 at midpoint as in Julia workflow for coupling matrix element reference
    Q0mid = float(poti.Q[len(poti.Q)//2])
    calc_overlap(cc, Q0=Q0mid)
    temps = np.linspace(Tmin, Tmax, NT)
    calc_capt_coeff(cc, V=V, temperature=temps)
    return cc


def get_morse_capture(a_i: float, a_f: float, b_i: float, b_f: float, Q0: float, E0: float, W: float, V: float, Tmin: float = 100.0, Tmax: float = 800.0, NT: int = 100):
    poti, potf = fit_morse_params(a_i, a_f, b_i, b_f, Q0, E0)
    cc = ConfCoord(V1=poti, V2=potf, W=W, g=1)
    Q0mid = float(poti.Q[len(poti.Q)//2])
    calc_overlap(cc, Q0=Q0mid)
    temps = np.linspace(Tmin, Tmax, NT)
    calc_capt_coeff(cc, V=V, temperature=temps)
    return cc
