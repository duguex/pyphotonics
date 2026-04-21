"""Plot helpers (ported from Plotter.jl) using matplotlib."""

from __future__ import annotations

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from .potential import Potential
from .capture_rate import ConfCoord


def plot_pot(pot: Potential, *, plot_wf: bool = False, ax: Optional[plt.Axes] = None, label: str = "", scale_factor: float = 2e-2):
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(pot.Q, pot.E, label=label or pot.name)

    if plot_wf and pot.chi.size:
        # offset wavefunctions by eigenvalues
        for i, eps in enumerate(pot.eps):
            wf = pot.chi[i, :]
            ax.plot(pot.Q, eps + scale_factor * wf, alpha=0.6)

    ax.set_xlabel("Q")
    ax.set_ylabel("Energy (eV)")
    return ax


def plot_cc(cc: ConfCoord, *, ax: Optional[plt.Axes] = None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(cc.temperature, cc.capt_coeff)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Capture coefficient (cm^3/s)")
    ax.set_yscale("log")
    return ax
