#!/share/home/ckduan/anaconda3/envs/my_pydefect/bin/python
from photonics2.photoluminescence import Photoluminescence

from photonics2.plott import plot_S_I
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def get_pl(path_phonopy, ground_struct_path, excited_struct_path, **parameter):
    res = parameter.get("resolution", 1000)
    proc = parameter.get("process", "emission")
    plot_width = parameter.get("plot_width", 1.0)
    if "emi" in proc:
        Amin = -0.8 * plot_width
        Amax = 0.2 * plot_width
    elif "abs" in proc:
        Amin = -0.2 * plot_width
        Amax = 0.8 * plot_width

    print("################# reading band.yaml Start! ###################")
    mass = parameter.get("mass", [])
    if len(mass) > 0:
        p = Photoluminescence(path_phonopy, "phonopy", m=mass, POSCAR_GRD=ground_struct_path,
                              POSCAR_EX=excited_struct_path,n_defect=1, resolution=res,
                              defect_range=-1.0)
    else:
        p = Photoluminescence(path_phonopy, "phonopy", POSCAR_GRD=ground_struct_path, POSCAR_EX=excited_struct_path,
                              n_defect=1, resolution=res, defect_range=-1.0)
    print("################# band.yaml OK! ###################")

    print("################# HuangRhyes Start! ###################")
    HR = p.HuangRhyes()
    print("HuangRhyes=", HR)

    f = open("partial.HuangRhyes.data", "w")
    f.write("mode No.\t partial.HuangRhyes \n")
    i = 0
    for s in p.S:
        i += 1
        f.write(str(i) + "\t" + str(s) + "\n")
    f.close()
    p.print_table(0.1, [])

    print("################# HuangRhyes OK! ###################")

    print("################# el_ph Start! ###################")
    gw = parameter.get("gw", 1e-3)

    p.el_ph(delta_width=gw, temperature=parameter.get("T", 0), jtmodes=[])
    print("################# el_ph OK! ###################")

    print("################# HuangRhyes lineshape Start! ###################")
    zpl = parameter.get("zpl", 2.5)

    A = p.PL(gamma=parameter.get("gamma", 10e-3), SHR=0, EZPL=zpl, process=proc)

    p.PLA()
    print("################# HuangRhyes lineshape OK! ###################")
    plot_S_I(p, parameter.get("title", "any"), [p.EZPL + Amin * 0.5, p.EZPL + Amax * 0.5], "Shw PLeV  PLnm ")

    return p


def read_parameters_from_incar(incar_path):
    params = {}
    with open(incar_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                try:
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
                params[key] = value
    return params


incar_path = ".\INCAR"
if len(sys.argv) > 1:
    incar_path = sys.argv[1]

incar_params = read_parameters_from_incar(incar_path)
get_pl(
    incar_params.get("BAND_YAML", "./zlq/band.yaml"),
    incar_params.get("GRD_POSCAR", "./zlq/GS"),
    incar_params.get("EXC_POSCAR", "./zlq/ES"),
    gamma=incar_params.get("gamma", 10e-3),
    zpl=incar_params.get("zpl", 1.0),
    gw=incar_params.get("gw", 1e-3),
    resolution=incar_params.get("resolution", 1000),
    T=incar_params.get("T", 0),
    plot_width=incar_params.get("plot_width", 2.5),
    title=incar_params.get("title", "any"),
    process=incar_params.get("process", "emission")
)