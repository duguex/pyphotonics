"""FERWE/FERDO INCAR tag computation for constrained DFT in VASP.

Provides the ``ferwe_ferdo()`` function which reads a VASP OUTCAR
and prints the INCAR tags needed to constrain electronic occupations.
"""

from __future__ import annotations

from typing import Tuple

from oganesson.io.vasp import Outcar


def ferwe_ferdo(vasp_folder: str = "./") -> Tuple[str, str]:
    """Read OUTCAR and return FERWE and FERDO INCAR tag strings.

    Parameters
    ----------
    vasp_folder : str
        Directory containing the VASP OUTCAR.

    Returns
    -------
    tuple of str
        ``(ferwe, ferdo)`` formatted INCAR tag values.
    """
    outcar = Outcar(vasp_folder, "OUTCAR")
    return outcar.get_ferwe_ferdo()
