"""CarrierCapture (Python port)

This package is a Python translation of the Julia repo in `src/`.
It provides:
  - Potential: 1D PES container + fitting + 1D Schrödinger solver
  - ConfCoord: configuration-coordinate capture coefficient utilities
  - TransferCoord: diabatic/adiabatic transfer rate utilities
  - Brooglie: sparse-grid Schrödinger solver backend (ported from Brooglie.jl)

The public API mirrors the Julia exports as closely as practical in Python.
"""

from .constants import HBAR_EV_S, KB_EV_K, AMU_EV_C2, HBARC_EV_M
from .brooglie import build_h, solve
from .potential import (
    Potential,
    pot_from_file,
    filter_sample_points,
    fit_pot,
    solve_pot,
    find_crossing,
    cleave_pot,
)
from .capture_rate import ConfCoord, calc_overlap, calc_capt_coeff
from .transfer_coord import TransferCoord, get_coupling, get_reorg_energy, get_activation_energy, get_transfer_rate, einstein_mobility

__all__ = [
    "HBAR_EV_S","KB_EV_K","AMU_EV_C2","HBARC_EV_M",
    "build_h","solve",
    "Potential","pot_from_file","filter_sample_points","fit_pot","solve_pot","find_crossing","cleave_pot",
    "ConfCoord","calc_overlap","calc_capt_coeff",
    "TransferCoord","get_coupling","get_reorg_energy","get_activation_energy","get_transfer_rate","einstein_mobility",
]
