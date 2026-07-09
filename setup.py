from setuptools import setup
from pathlib import Path

here = Path(__file__).resolve().parent
README = (here / "README.md").read_text(encoding="utf-8")
VERSION = (here / "pyphotonics" / "VERSION").read_text(encoding="utf-8").strip()

setup(
    name="pyphotonics",
    packages=["pyphotonics", "carriercapture"],
    entry_points={
        "console_scripts": [
            "pyphotonics=pyphotonics.cli:execute_cli",
            "pyphotonics-incar=pyphotonics.cli:execute_incars",
        ],
    },
    include_package_data=True,
    version=VERSION,
    license="gpl-3.0",
    description="Post-processing Python code that calculates photonic properties of defects using output files from VASP and phonopy. Computes Huang-Rhys factor, photoluminescence line-shapes, and carrier capture coefficients.",
    author="Sherif Abdulkader Tawfik",
    author_email="sherif.tawfic@gmail.com",
    long_description=README,
    long_description_content_type='text/markdown',
    url="https://github.com/sheriftawfikabbas/pyphotonics",
    keywords=["DFT", "Material science", "Photoluminescence", "VASP", "Carrier capture"],
    install_requires=[
        "scipy",
        "numpy",
        "pandas",
        "matplotlib",
        "pymatgen",
        "oganesson",
    ],
)
