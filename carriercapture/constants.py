"""Physical constants used throughout the port (same numeric values as the Julia code)."""

import math

# Potential.jl constants
AMU_EV_C2 = 931.4940954e6   # eV / c^2
HBARC_EV_M = 0.19732697e-6  # eV·m
BOLTZ_EV_K = 8.617333262e-5 # eV/K

# CaptureRate.jl / TransferCoord.jl constants
HBAR_EV_S = 6.582119514e-16  # eV·s
KB_EV_K = 8.6173303e-5       # eV·K^-1 (kept as in Julia; close to BOLTZ_EV_K)

EV_TO_JOULE = 1.60218e-19
E_CHARGE_C = 1.60217662e-19

TWO_PI = 2.0 * math.pi
