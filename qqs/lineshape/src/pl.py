#!/home/qqs/miniconda3/envs/pymatgen/bin/python
from photonics2.photoluminescence import Photoluminescence

from photonics2.plott import plot_S_I
import numpy as np
import os
import sys
#import photonics.jt as jt
import matplotlib.pyplot as plt

def get_pl(path_phonopy, ground_struct_path, excited_struct_path, **parameter):
    res = 1000
    proc= parameter.get("process", "emission")  # emission or absorption
    plot_width= parameter.get("plot_width", 1.0)  # width of the plot
    if "emi" in proc:
        Amin=-0.9*plot_width
        Amax=0.1*plot_width
    elif "abs" in proc:
        Amin=-0.1*plot_width
        Amax=0.9*plot_width
    ############ Jahn-Teller ##############

    ################# HR lumine ###################
    print("################# reading band.yaml Start! ###################")
    mass = parameter.get("mass", [])
    if len(mass) > 0:
        p = Photoluminescence(path_phonopy, "phonopy", m=mass, POSCAR_GRD=ground_struct_path,
                              POSCAR_EX=excited_struct_path, n_defect=1, resolution=res)
    else:
        p = Photoluminescence(path_phonopy, "phonopy", POSCAR_GRD=ground_struct_path, POSCAR_EX=excited_struct_path,
                              n_defect=1, resolution=res)
    print("################# band.yaml OK! ###################")
    # p = Photoluminescence(path_phonopy,"phonopy large",m=masscsgo,POSCAR_GRD="./CsSO-pbesol/HS",\
    # force_grd="./CsSO-pbesol/hs_force.test",force_ex="./CsSO-pbesol/ls2_force.test"\
    # ,n_defect=1,resolution=res,  shift_vector=[0.0, 0.0, 0.0])

    # Photoluminescence can accept use's input of mass list (m=[m1,m2....] and  nomber of  modes (nmodes = int),then it reads eigenvector and frequences automatically.
    # or Photoluminescence accepts none of these imformation, it will read all these needed data automatically but inefficiently.
    # p = Photoluminescence(path_phonopy,"phonopy",POSCAR_GRD=ground_struct_path,POSCAR_EX=excited_struct_path,n_defect=1,resolution=res,  shift_vector=[0.1, 0.1, 0.0])
    # print(p.m)
    
    

    ################# HuangRhyes ###################

    print("################# HuangRhyes Start! ###################")
    
    #p.frequencies[384]=0.04039
    #p.frequencies[385]=0.04021
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
    ############# #el_ph coupling ################ 
    print("################# el_ph Start! ###################")
    gw =parameter.get("gw", 1e-3)

    p.el_ph(delta_width=gw, temperature=parameter.get("T", 0), jtmodes=[])
    print("################# el_ph OK! ###################")

    
    ################# HuangRhyes lineshape################
    print("################# HuangRhyes lineshape Start! ###################")
    zpl = parameter.get("zpl", 2.5)
    #zpl = 1.8391+0.0439+0.038+0.033+(720-1900)*0.000124-0.0855-0.208-0.01+0.0292#abs
    #zpl = 10557 * 0.000124
    
    #print( 0.8391+0.0429-10557*0.000124)
    A = p.PL(gamma=parameter.get("gamma", 10e-3)*res, SHR=0, EZPL=zpl, process=proc)  # gamma ZHR zpl  ZHR only influences compute
    #plot_S_I(p, "$K_2SO_4:Mn^{6+}$", [6200 * 0.000124, 11100 * 0.000124], "Shw+Acm PLcm",split=122.817*0.000124)
    
    # input:gamma(for width of zpl) SHR(for calculation) EZPL(for energy of zpl) p rocess(emission or absorption)
    
    ############## convolve all line ##################
    #plot_S_I(p, "$Cs_2SO_4:Mn^{6+} $ ", [8400 * 0.000124, 10700 * 0.000124], " Shw+Acm ")

    p.PLA()
    
    print("################# HuangRhyes lineshape OK! ###################")
    plot_S_I(p, parameter.get("title", "any"), [p.EZPL+Amin*0.5 , p.EZPL+Amax*0.5], "Shw Sk PLeV  PLnm A Phon")
    #plot_S_I(p, "$^2T_2 $ ", [10300 * 0.000124, 10700 * 0.000124], "  PLcm ",)
    
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
                # Try to convert to float or int if possible
                try:
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
                params[key] = value
    return params

# 默认参数文件路径
incar_path = os.path.join(os.path.dirname(__file__), "INCAR")
# 如果命令行有参数，则使用指定路径
if len(sys.argv) > 1:
    incar_path = sys.argv[1]

incar_params = read_parameters_from_incar(incar_path)
get_pl(
    incar_params.get("BAND_YAML", "../cases/zlq/band.yaml"),
    incar_params.get("GRD_POSCAR", "../cases/zlq/GS"),
    incar_params.get("EXC_POSCAR", "../cases/zlq/ES"),
    gamma=incar_params.get("gamma", 10e-3),
    zpl=incar_params.get("zpl", 1.0),
    gw=incar_params.get("gw", 1e-3),
    resolution=incar_params.get("resolution", 1000),
    T=incar_params.get("T", 0),
    plot_width=incar_params.get("plot_width", 2.5),
    title=incar_params.get("title", "any"),
    process=incar_params.get("process", "emission")
)



   
"""
Calculates and plots the photoluminescence (PL) spectrum using the Photoluminescence class from the photonics2 library.

This function reads phonon and structural data, computes Huang-Rhys factors, electron-phonon coupling, and generates the PL lineshape.
It also saves partial Huang-Rhys factors to a file and produces plots of the resulting spectra.

Parameters:
    path_phonopy (str): Path to the phonopy band.yaml file.
    ground_struct_path (str): Path to the ground state structure (POSCAR file).
    excited_struct_path (str): Path to the excited state structure (POSCAR file).
    **parameter: Additional keyword arguments:
        - process (str): 'emission' or 'absorption' (default: 'emission').
        - gamma (float): Broadening parameter for the ZPL (default: 10e-3). (eV)
        - zpl (float): Zero-phonon line energy (default: 2.5). (eV)
        - gw (float): Width parameter for electron-phonon coupling (default: 1e-3). (eV)
        - resolution (int): Number of points in the spectrum (default: 1000). (1000 points corresponds to 1 eV.)
        - plot_width (float): Width of the plot window (default: 1.0). （eV）
        - T (float): Temperature for electron-phonon coupling (default: 0).
        - title (str): Title for the plot (default: "any").

Returns:
    Photoluminescence: An instance of the Photoluminescence class with computed properties.

Side Effects:
    - Writes partial Huang-Rhys factors to "partial.HuangRhyes.data".
    - Prints progress and results to the console.
    - Generates and displays plots of the PL spectrum. Shw.png Sk.png PLeV.png  PLnm.png A.png Phon.png.    

Example:
    get_pl(
        "./Cs3Cu2Br5_STE/band.yaml",
        "./Cs3Cu2Br5_STE/GS",
        "./Cs3Cu2Br5_STE/ES",
        gamma=10e-3, zpl=3.24, gw=1e-3, resolution=1000, T=0, plot_width=2.0, process="emission"
    )
"""
