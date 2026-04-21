from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import HBAR_EV_S, KB_EV_K, EV_TO_JOULE, E_CHARGE_C
from .potential import Potential, find_crossing


@dataclass
class TransferCoord:
    """Three potentials describing a transfer: two diabatic and one adiabatic."""
    pot_d_i: Potential
    pot_d_f: Potential
    pot_a: Potential


def get_coupling(tc: TransferCoord) -> float:
    crossing_Q, crossing_E = find_crossing(tc.pot_d_i, tc.pot_d_f)
    coupling = crossing_E - float(tc.pot_a.func([crossing_Q])[0])
    return float(coupling)


def get_reorg_energy(tc: TransferCoord) -> float:
    Q_min_i = tc.pot_d_i.Q0
    E_min_f = tc.pot_d_f.E0
    E_min_i_f = float(tc.pot_d_f.func([Q_min_i])[0])
    return float(E_min_i_f - E_min_f)


def get_activation_energy(tc: TransferCoord) -> float:
    _Qx, Ex = find_crossing(tc.pot_d_i, tc.pot_d_f)
    return float(Ex - tc.pot_d_i.E0)


def get_transfer_rate(tc: TransferCoord) -> float:
    lam = get_reorg_energy(tc)
    Hab = get_coupling(tc)
    dG = get_activation_energy(tc)
    kBT = KB_EV_K * tc.pot_a.T

    rate = (2.0 * math.pi / HBAR_EV_S) * (1.0 / math.sqrt(4.0 * math.pi * lam * kBT)) * (Hab ** 2) * math.exp(-dG / kBT)
    return float(rate)


def einstein_mobility(rate: float, n_neighbours: int, dist: float, temp: float) -> float:
    diffusion = (dist ** 2) * int(n_neighbours) * float(rate)
    mob = E_CHARGE_C * diffusion / (KB_EV_K * EV_TO_JOULE * float(temp))
    return float(mob)
