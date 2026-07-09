"""PyPhotonics — post-processing DFT code for photonic properties of defects.

Subpackages:
  pyphotonics          — Photoluminescence line-shapes (Huang-Rhys factor, PL spectra)
  carriercapture       — Carrier capture coefficients, 1D PES, Schrödinger solver
"""

from .version import VERSION
from .photoluminescence import Photoluminescence

__all__ = [
    "VERSION",
    "Photoluminescence",
]
