"""Command-line interfaces for PyPhotonics."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

from oganesson.io.vasp import Outcar

from .photoluminescence import Photoluminescence
from .version import VERSION


class CLI_Photoluminescence:
    """CLI for Huang-Rhys factor and PL line-shape calculation."""

    def __init__(self, argv: Optional[list[str]] = None) -> None:
        argv = argv or sys.argv[:]
        prog = Path(argv[0]).name

        parser = argparse.ArgumentParser(
            prog=prog,
            description=f"{prog} version {VERSION}",
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument(
            "--version", action="version", version=f"%(prog)s {VERSION}"
        )
        parser.add_argument(
            "-cgs",
            "--contcar_ground_state",
            type=str,
            required=True,
            help="filename of the ground-state CONTCAR",
        )
        parser.add_argument(
            "-ces",
            "--contcar_excited_state",
            type=str,
            required=True,
            help="filename of the excited-state CONTCAR",
        )
        parser.add_argument(
            "-m",
            "--num_modes",
            type=int,
            required=True,
            help="number of vibrational modes",
        )
        parser.add_argument(
            "-M",
            "--method",
            type=str,
            default="phonopy",
            help='method: "phonopy", "phonopy-siesta", or "vasp"',
        )
        parser.add_argument(
            "-r",
            "--resolution",
            type=int,
            default=1000,
            help="energy grid resolution (points per eV)",
        )
        parser.add_argument(
            "--phonopy_path",
            type=str,
            default="./phonopy",
            help="directory containing phonopy band.yaml",
        )

        args = parser.parse_args(argv[1:])

        pl = Photoluminescence(
            ground_state=args.contcar_ground_state,
            excited_state=args.contcar_excited_state,
            num_modes=args.num_modes,
            method=args.method,
            resolution=args.resolution,
            phonopy_path=args.phonopy_path,
        )
        pl.run()


class CLI_INCARs:
    """CLI for generating FERWE/FERDO INCAR tags for constrained DFT."""

    def __init__(self, argv: Optional[list[str]] = None) -> None:
        argv = argv or sys.argv[:]
        prog = Path(argv[0]).name

        parser = argparse.ArgumentParser(
            prog=prog,
            description=f"{prog} version {VERSION}",
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument(
            "--version", action="version", version=f"%(prog)s {VERSION}"
        )
        parser.add_argument(
            "-f",
            "--vasp_folder",
            type=str,
            default="./",
            help="directory containing the VASP OUTCAR file",
        )

        args = parser.parse_args(argv[1:])

        outcar = Outcar(args.vasp_folder, "OUTCAR")
        up, down = outcar.get_ferwe_ferdo()
        print("FERWE =", up)
        print("FERDO =", down)


def execute_cli(argv: Optional[list[str]] = None) -> None:
    CLI_Photoluminescence(argv)


def execute_incars(argv: Optional[list[str]] = None) -> None:
    CLI_INCARs(argv)
