"""Photoluminescence line-shape calculation for defects in crystals.

Computes the Huang-Rhys factor and PL line-shape from VASP + phonopy output.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, Union

import numpy as np
from numpy import fft
from oganesson import OgStructure
from pymatgen.core import Structure

# ── Physical constants ────────────────────────────────────────────────
HBAR_JS = 1.054571817e-34          # ħ in J·s
HBAR_EVS = 6.582119569e-16         # ħ in eV·s
AMU_TO_KG = 1.66053906660e-27      # amu → kg
THZ_TO_EV = 0.004135665538536       # THz → eV
SIESTA_FACTOR = 0.727445665         # phonopy-siesta correction
DEFAULT_MAX_ENERGY = 5.0            # eV
DEFAULT_SIGMA = 6e-3                # eV Gaussian broadening
DEFAULT_RESOLUTION = 1000


# ── Public class ──────────────────────────────────────────────────────

class Photoluminescence:
    """Compute Huang-Rhys factor and PL line-shape for a defect system.

    Reads ground/excited structures (VASP CONTCAR), phonon modes
    (phonopy ``band.yaml``), and computes the Huang-Rhys factor per
    vibrational mode, the spectral function S(ω), and the PL line-shape.

    Parameters
    ----------
    ground_state : str
        Ground-state CONTCAR filename (relative to ``path``).
    excited_state : str
        Excited-state CONTCAR filename (relative to ``path``).
    num_modes : int
        Number of vibrational modes.
    method : str
        ``"phonopy"``, ``"phonopy-siesta"``, or ``"vasp"``.
    path : str
        Base directory for CONTCAR files.
    phonopy_path : str
        Directory containing ``band.yaml``.
    masses : np.ndarray | None
        Atomic masses in kg (auto-computed from structure if ``None``).
    resolution : int
        Energy grid resolution (points per eV).

    Attributes
    ----------
    HuangRhys : float
        Total Huang-Rhys factor.
    Delta_R : float
        RMS structural displacement (Å).
    Delta_Q : float
        Mass-weighted displacement (amu^1/2 . Angstrom).
    frequencies : np.ndarray
        Mode frequencies (eV), shape ``(num_modes,)``.
    S : list[float]
        Mass-weighted displacement (amu^1/2 · Å).
    q : list[float]
        Displacement coordinate per mode.
    IPR : list[float]
        Inverse participation ratio per mode.
    S_omega : list[float]
        Spectral Huang-Rhys function on energy grid.
    omega_set : np.ndarray
        Energy grid for ``S_omega`` (eV).
    """

    # ── Construction ──────────────────────────────────────────────

    def __init__(
        self,
        ground_state: str,
        excited_state: Optional[str] = None,
        num_modes: Optional[int] = None,
        method: str = "phonopy",
        path: str = "./",
        phonopy_path: str = "./phonopy",
        masses: Optional[np.ndarray] = None,
        resolution: int = DEFAULT_RESOLUTION,
        # ── backward-compat aliases ──
        exceited_state: Optional[str] = None,
        numModes: Optional[int] = None,
        m: Optional[np.ndarray] = None,
    ) -> None:
        # Resolve backward-compat names
        _excited_state = excited_state or exceited_state
        _num_modes = num_modes if num_modes is not None else numModes
        _masses = masses if masses is not None else m

        if _excited_state is None:
            raise TypeError("missing required argument: excited_state")
        if _num_modes is None:
            raise TypeError("missing required argument: num_modes")

        self.path = Path(path)
        self.phonopy_path = Path(phonopy_path)
        self.resolution = resolution
        self.num_modes = _num_modes
        self.method = method

        # ── Load structures ──────────────────────────────────────
        gs_file = self.path / ground_state
        es_file = self.path / _excited_state

        if not gs_file.exists():
            raise FileNotFoundError(f"ground state not found: {gs_file}")
        if not es_file.exists():
            raise FileNotFoundError(f"excited state not found: {es_file}")

        # Read both structures with the same primitive (pymatgen), so
        # the per-atom index in `delta_r` aligns with the Modes array
        # (which is read from band.yaml in POSCAR input order). Using
        # OgStructure here would re-order atoms by electronegativity
        # and break the index alignment — see tools/HR_DIVERGENCE.md.
        self._ground_structure = Structure.from_file(str(gs_file))
        self._excited_structure = Structure.from_file(str(es_file))
        from pymatgen.core.lattice import pbc_shortest_vectors
        delta_r = pbc_shortest_vectors(
            self._ground_structure.lattice,
            self._ground_structure.frac_coords,
            self._excited_structure.frac_coords,
        )[np.arange(len(self._ground_structure)),
          np.arange(len(self._excited_structure))]
        self.num_atoms = len(self._ground_structure)

        # ── Masses ───────────────────────────────────────────────
        if _masses is not None:
            self._masses = np.asarray(_masses, dtype=float)
        else:
            self._masses = np.array([
                s.atomic_mass * AMU_TO_KG
                for s in self._ground_structure.species
            ], dtype=float)

        # ── Load phonon data ─────────────────────────────────────
        if "phonopy" in method:
            self.frequencies = self._phonopy_read_frequencies()
            modes = self._phonopy_read_modes()
        else:
            raise NotImplementedError(
                f"method={method!r} is not implemented; use 'phonopy'"
            )

        # ── Compute ──────────────────────────────────────────────
        self._compute_hr(modes, delta_r)
        self._compute_s_omega()

    # ── Phonon file readers ──────────────────────────────────────

    def _phonopy_read_modes(self) -> np.ndarray:
        """Read eigen-displacement vectors from ``band.yaml``.

        Returns
        -------
        np.ndarray
            Shape ``(num_modes, num_atoms, 3)``.
        """
        modes = np.zeros((self.num_modes, self.num_atoms, 3))
        path = self.phonopy_path / "band.yaml"

        try:
            with open(path) as fh:
                for line in fh:
                    if "  band:" in line:
                        break

                for i in range(self.num_modes):
                    fh.readline()  # - band #
                    fh.readline()  # frequency: ...
                    fh.readline()  # eigenvector:
                    for a in range(self.num_atoms):
                        fh.readline()  # - # atom N  (consume label)

                        line = fh.readline().replace(",", "")
                        parts = line.strip().split()
                        modes[i][a][0] = float(parts[2])

                        line = fh.readline().replace(",", "")
                        parts = line.strip().split()
                        modes[i][a][1] = float(parts[2])

                        line = fh.readline().replace(",", "")
                        parts = line.strip().split()
                        modes[i][a][2] = float(parts[2])
        except OSError as exc:
            raise FileNotFoundError(f"could not read {path}") from exc

        return modes

    def _phonopy_read_frequencies(self) -> np.ndarray:
        """Read phonon frequencies from ``band.yaml``.

        Returns
        -------
        np.ndarray
            Shape ``(num_modes,)``, frequencies in eV.
        """
        freqs = np.zeros(self.num_modes)
        path = self.phonopy_path / "band.yaml"

        try:
            with open(path) as f:
                for line in f:
                    if "  band:" in line:
                        break

                for i in range(self.num_modes):
                    f.readline()  # - band #:
                    line = f.readline()
                    parts = line.strip().split()
                    freqs[i] = float(parts[1])
                    f.readline()  # atom:
                    for _ in range(self.num_atoms):
                        f.readline()  # displacement:
                        f.readline()  # x
                        f.readline()  # y
                        f.readline()  # z
        except OSError as exc:
            raise FileNotFoundError(f"could not read {path}") from exc

        # Convert units and clip
        if self.method == "phonopy":
            freqs = freqs * THZ_TO_EV
        elif self.method == "phonopy-siesta":
            freqs = freqs * THZ_TO_EV * SIESTA_FACTOR
        # NOTE: "vasp" method not implemented

        freqs[freqs < 0] = 0.0
        return freqs

    # ── Core computation ─────────────────────────────────────────

    def _compute_hr(self, modes: np.ndarray, delta_r: np.ndarray) -> None:
        """Compute Huang-Rhys factor, Δ_R, Δ_Q, IPR, q per mode.

        This reproduces the original computation exactly.
        """
        num_modes = self.num_modes
        num_atoms = self.num_atoms

        self.HuangRhys = 0.0
        self.Delta_R = 0.0
        self.Delta_Q = 0.0
        self.IPR: List[float] = []
        self.q: List[float] = []
        self.S: List[float] = []

        for i in range(num_modes):
            q_i = 0.0
            ipr_i = 0.0

            for a in range(num_atoms):
                rx, ry, rz = modes[i][a]
                participation = rx * rx + ry * ry + rz * rz
                ipr_i += participation ** 2

                q_i += (
                    np.sqrt(self._masses[a])
                    * delta_r[a][0] * rx
                    + np.sqrt(self._masses[a])
                    * delta_r[a][1] * ry
                    + np.sqrt(self._masses[a])
                    * delta_r[a][2] * rz
                ) * 1e-10

            ipr_i = 1.0 / ipr_i if ipr_i > 0 else float("inf")

            # Huang-Rhys factor for this mode
            # ω [eV] · q² [kg·m²] / 2 / (ħ [J·s] · ħ [eV·s])
            omega = float(self.frequencies[i])
            s_i = omega * q_i * q_i / 2.0 / (HBAR_JS * HBAR_EVS)

            self.IPR.append(ipr_i)
            self.q.append(q_i)
            self.S.append(s_i)
            self.HuangRhys += s_i

        # Δ_R and Δ_Q
        for a in range(num_atoms):
            for coord in range(3):
                d = delta_r[a][coord]
                self.Delta_R += d * d
                self.Delta_Q += d * d * self._masses[a]

        self.Delta_R = float(np.sqrt(self.Delta_R))
        self.Delta_Q = float(np.sqrt(self.Delta_Q / AMU_TO_KG))

    def _compute_s_omega(self) -> None:
        """Compute spectral function S(ω) on an energy grid."""
        max_energy = DEFAULT_MAX_ENERGY
        sigma = DEFAULT_SIGMA

        self.omega_set = np.linspace(
            0, max_energy, int(max_energy * self.resolution)
        )
        self.S_omega = [
            self.get_S_omega(omega, sigma) for omega in self.omega_set
        ]

    # ── Public API ───────────────────────────────────────────────

    def get_S_omega(self, omega: float, sigma: float) -> float:
        """Spectral Huang-Rhys function at a single energy.

        Parameters
        ----------
        omega : float
            Energy (eV).
        sigma : float
            Gaussian broadening width (eV).

        Returns
        -------
        float
            S(ω) value.
        """
        total = 0.0
        for s_k, freq in zip(self.S, self.frequencies):
            total += s_k * self._gaussian(omega, freq, sigma)
        return total

    @staticmethod
    def _gaussian(omega: float, omega_k: float, sigma: float) -> float:
        """Normalised Gaussian centred on ``omega_k``."""
        arg = (omega - omega_k) / sigma
        return np.exp(-0.5 * arg * arg) / (np.sqrt(2.0 * np.pi) * sigma)

    def write_S(self, file_name: Union[str, Path]) -> None:
        """Write S(ω) to a text file, one value per line."""
        with open(file_name, "w") as f:
            for val in self.S_omega:
                f.write(str(val) + "\n")

    def PL(
        self,
        gamma: float,
        SHR: float,
        EZPL: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute photoluminescence line-shape.

        Uses Fourier transform of the spectral function.

        Parameters
        ----------
        gamma : float
            Lifetime broadening (eV).
        SHR : float
            Total Huang-Rhys factor.
        EZPL : float
            Zero-phonon line energy (eV).

        Returns
        -------
        A : np.ndarray
            Absorption/emission kernel.
        I : np.ndarray
            PL intensity (``A × ω³``).
        """
        r = 1.0 / self.resolution

        St = fft.ifft(np.asarray(self.S_omega, dtype=complex))
        St = fft.ifftshift(St)
        G = np.exp(2.0 * np.pi * St - float(SHR))

        n = len(G)
        # Apply gamma broadening as a Lorentzian convolution in the
        # energy domain instead of exp(-gamma*|t|) in the time domain.
        # This decouples the line-shape width from `resolution` —
        # the convolution kernel depends only on `gamma`, not on the
        # FFT grid spacing. See tools/RESOLUTION_NOTES.md for analysis.
        A = fft.fft(G)
        omega = np.arange(n) * r - (n // 2) * r
        kernel = (gamma / np.pi) / (omega**2 + gamma**2)
        kernel = kernel / kernel.sum()
        A = np.convolve(A, kernel, mode="same")

        # Shift ZPL peak to EZPL
        tA = A.copy()
        shift = int(EZPL * self.resolution)
        for i in range(n):
            A[(shift - i) % n] = tA[i]

        I = np.array([A[i] * (i * r) ** 3 for i in range(n)])
        return A, I

    def print_table(self) -> None:
        """Print IPR, S_k, and frequency per mode to stdout."""
        for i in range(self.num_modes):
            print(
                f"IPR\t{i}\tSk\t{self.S[i]}\tenergy\t"
                f"{self.frequencies[i]}\t=\t{self.IPR[i]}"
                f"\twith localization ratio beta =\t{64 / self.IPR[i]}"
            )

    def run(self) -> None:
        """Entry point called by the CLI."""
        return  # all work done in __init__

    # ── Backward-compat accessors ───────────────────────────────
    @property
    def HuangRhyes(self) -> float:
        """Backward-compat alias for ``HuangRhys`` (original typo)."""
        return self.HuangRhys

    @property
    def numAtoms(self) -> int:
        """Backward-compat alias for ``num_atoms``."""
        return self.num_atoms

