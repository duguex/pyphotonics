import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator


def plot_S_I(p, name, plrange, image, **para):
    print(image)

    if "q" in image:
        fig = plt.figure(figsize=(8, 4), dpi=500)
        ax = fig.add_subplot()
        ax.set_title(name + " (total Q = " + str(round(p.Delta_Q, 3)) + "amu^1/2*A)")
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_xlim(0, max(p.frequencies) * 1.1)
        markerline, stemlines, baseline = ax.stem(p.frequencies[p.skipmodes + 1:], p.q_amu[p.skipmodes + 1:], label="partial dQ", basefmt="k-")
        markerline.set_markerfacecolor("tan")
        ax.set_ylabel("partial DQ")
        ax.legend(loc="upper right")
        ax.set_ylim(min(p.q_amu[p.skipmodes + 1:]) * 1.2, max(p.q_amu[p.skipmodes + 1:]) * 1.2)
        plt.savefig("q.png", dpi=500)

    if "Sk" in image:
        fig = plt.figure(figsize=(8, 4), dpi=500)
        ax = fig.add_subplot()
        ax.set_title(name + " (total HR = " + str(round(p.HuangRhyes, 3)) + ")")
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_xlim(0, max(p.frequencies) * 1.1)
        markerline, stemlines, baseline = ax.stem(p.frequencies[p.skipmodes + 1:], p.S[p.skipmodes + 1:], label="partial HuangRhyes", basefmt="k-")
        markerline.set_markerfacecolor("tan")
        ax.set_ylabel("partial HuangRhyes Sk")
        ax.legend(loc="upper right")
        ax.set_ylim(0, max(p.S[p.skipmodes + 1:]) * 1.2)
        plt.savefig("Sk.png", dpi=500)

    if "Shw" in image:
        plt.rcParams['font.size'] = 12
        fig = plt.figure(figsize=(6, 5), dpi=300)
        ax = fig.add_subplot()
        plt.subplots_adjust(bottom=0.15, top=0.95)
        ax2 = ax.twinx()
        markerline, stemlines, baseline = ax2.stem(p.frequencies[p.skipmodes:], p.S[p.skipmodes:], label="$S_k$", basefmt="k-")
        ax2.set_ylabel("$S_k$")
        ax2.set_ylim(0, max(p.S) * 1.2)

        ax.plot(p.omega_set, p.C_omega, color="tomato", label="C($\hbar\omega,T$) ")
        ax.plot(p.omega_set, p.S_omega, color="aqua", label="S($\hbar\omega$)")
        data_to_save = np.column_stack((p.omega_set, p.C_omega, p.S_omega))
        np.savetxt('C_omega,S_omega.data', data_to_save, delimiter=',', header='omega,C_omega,S_omega', comments='')
        ax.set_ylabel('S($\hbar\omega$) (1/eV)')
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_xlim(0, max(p.frequencies) * 1.1)
        ax.set_ylim(0, max([max(p.S_omega), max(p.C_omega)]) * 1.5)
        fig.legend(loc=(0.70, 0.73))
        plt.savefig("Shw.png", dpi=500)

    if "Phon" in image:
        fig = plt.figure(figsize=(8, 4), dpi=500)
        ax = fig.add_subplot()
        ax.plot(p.omega_set, p.S_omega, color="teal", label="S($\hbar\omega$)  ")
        ax.set_title('Phonon DOS')
        ax.set_ylabel('DOS($\hbar\omega$) (1/eV)')
        ax.set_xlabel('Phonon energy (eV)')
        ax.set_xlim(0, max(p.frequencies) * 1.1)
        ax.set_ylim(0, max([max(p.S_omega), max(p.C_omega)]) * 1.2)
        ax.legend(loc="upper left")
        plt.savefig("Phon.png", dpi=500)

    if "A" in image:
        plt.figure(figsize=(8, 4), dpi=500)
        x = np.linspace(0, float(len(p.A)) / float(p.resolution), num=len(p.A))
        Aabs = p.A.__abs__() / sum(p.A.__abs__()) * p.resolution
        if p.process == "absorption":
            plt.plot(x, Aabs, label="PLE cross section", color="darkolivegreen")
            plt.title(name + "  PLE cross section A")
            plt.xlim(plrange[0], plrange[1])
            plt.ylabel('Cross section($\hbar\omega$)')
        elif p.process == "emission":
            plt.plot(x, Aabs, label="$A_e(\hbar\omega)$", color="darkolivegreen")
            plt.title(name + "  PL spectral function $A_e(\hbar\omega)$ ")
            plt.xlim(plrange[0], plrange[1])
            plt.ylabel('$Intensity(\hbar\omega)$')

        plt.xlabel('energy (eV)')
        plt.ylim(0, max(Aabs) * 1.1)
        plt.vlines(p.EZPL, ymin=0, ymax=max(Aabs) * 1.1, label='$E_{zpl}$=' + str(round(p.EZPL)) + "$eV$", color='grey', linestyle="--")
        plt.legend()
        plt.savefig("AeV.png", dpi=500)
        f = open("AeV.data", "w")
        i = 0
        for s in Aabs:
            f.write(str(x[i]) + "\t" + str(s) + "\n")
            i += 1
        f.close()

    if "Acm" in image:
        plt.figure(figsize=(8, 4), dpi=500)
        x = np.linspace(0, float(len(p.A)) / float(p.resolution), num=len(p.A))
        xcm = x / 0.000124
        Aabs = p.A.__abs__() / sum(p.A.__abs__()) * p.resolution * 0.000124
        if p.process == "absorption":
            plt.plot(xcm, Aabs, label="PLE cross section", color="darkolivegreen")
            plt.title(name + "  PLE cross section A")
            plt.xlim(plrange[0] / 0.000124, plrange[1] / 0.000124)
            plt.ylabel('Cross section($\hbar\omega$)')
        elif p.process == "emission":
            plt.plot(xcm, Aabs, label="$A_e(\hbar\omega)$", color="darkolivegreen")
            plt.title(name + "  PL spectral function $A_e(\hbar\omega)$ ")
            plt.xlim(plrange[0] / 0.000124, plrange[1] / 0.000124)
            plt.ylabel('Intensity$(\hbar\omega)$')

        plt.xlabel('energy (cm$^{-1}$)')
        plt.ylim(0, max(Aabs) * 1.1)
        plt.vlines(p.EZPL / 0.000124, ymin=0, ymax=max(Aabs) * 1.1, label='$E_{zpl}$=' + str(round(p.EZPL / 0.000124)) + "$cm^{-1}$", color='grey', linestyle="--")
        plt.legend()
        plt.savefig("Acm.png", dpi=500)

    if "Shw+Acm" in image:
        fig = plt.figure(figsize=(7, 7), dpi=500)
        ax = fig.add_subplot(2, 1, 1)
        ax.plot(p.omega_set, p.S_omega, linewidth=1.5, color="teal", label="S($\hbar\omega$)  ")
        ax.set_ylabel('S($\hbar\omega$) (1/eV)')
        ax.set_title(name)
        ax.set_xlim(0, max(p.frequencies) * 1.1)
        ax.set_ylim(0, max([max(p.S_omega), max(p.C_omega)]) * 1.2)
        ax.legend(loc="upper left")
        ax2 = ax.twinx()
        markerline, stemlines, baseline = ax2.stem(p.frequencies[p.skipmodes:], p.S[p.skipmodes:], label="partial HuangRhyes", basefmt="k-")
        markerline.set_markerfacecolor("tan")
        ax2.set_ylabel("partial HuangRhyes Sk")
        ax2.legend(loc="upper right")
        ax2.set_ylim(0, max(p.S) * 1.2)
        plt.savefig("Shw.png", dpi=500)

        ax3 = fig.add_subplot(2, 1, 2)
        x = np.linspace(0, float(len(p.A)) / float(p.resolution), num=len(p.A))
        xcm = x / 0.000124
        Aabs = p.A.__abs__() / sum(p.A.__abs__()) * p.resolution
        if p.process == "absorption":
            ax3.plot(x, Aabs, linewidth=1.5, label="PLE cross section", color="darkolivegreen")
            ax3.set_xlim(plrange[0], plrange[1])
            ax3.set_ylabel('Cross section($\hbar\omega$)')
        elif p.process == "emission":
            ax3.plot(x, Aabs, linewidth=1.5, label="$A_{other}(\hbar\omega)$", color="darkolivegreen")
            ax3.set_xlim(plrange[0], plrange[1])
            ax3.set_ylabel('Intensity$(\hbar\omega)$')

        ax3.set_xlabel('energy (eV)')
        ax3.set_ylim(0, max(Aabs) * 1.1)
        ax3.vlines(p.EZPL, ymin=0, ymax=max(Aabs) * 1.1, linewidth=1.5, label='$E_{zpl}$=' + str(round(p.EZPL, 3)) + "$eV$", color='grey', linestyle="--")
        ax3.legend()
        fig.savefig("Shw+Acm.png", dpi=500)
        fig.show()

    if "PLeV" in image:
        plt.figure(figsize=(8, 4), dpi=500)
        x = np.linspace(0, float(len(p.I)) / float(p.resolution), num=len(p.I))
        Iabs = p.I.__abs__() / sum(p.I.__abs__()) * p.resolution
        if p.process == "absorption":
            plt.plot(x, Iabs, label="PLE cross section", color="darkolivegreen")
            plt.title(name + "  PLE cross section")
            plt.xlim(plrange[0], plrange[1])
            plt.ylabel('Cross section($\hbar\omega$)')
        elif p.process == "emission":
            plt.plot(x, Iabs, label=" PL intensity", color="darkolivegreen")
            plt.title(name + "   PL intensity")
            plt.xlim(plrange[0], plrange[1])
            plt.ylabel(' PL Intensity')
        data_ple = np.column_stack((x, Iabs, Iabs / max(Iabs)))
        np.savetxt('lineshape_eV.csv', data_ple, delimiter=',', header='energy(eV),Iabs(a.u.),Iabs_norm(max=1)', comments='')
        plt.xlabel('energy (eV) T=' + str(p.temperature) + "K")
        plt.ylim(0, max(Iabs) * 1.1)
        plt.vlines(p.EZPL, ymin=0, ymax=max(Iabs) * 1.1, label='$E_{zpl}$=' + str(round(p.EZPL, 3)) + "eV", color='grey', linestyle="--")
        plt.legend()
        plt.savefig("PLev.png", dpi=500)
        f = open("PLev.data", "w")
        i = 0
        for s in Iabs:
            f.write(str(x[i]) + "\t" + str(s) + "\n")
            i += 1
        f.close()

    if "PLcm" in image:
        plt.figure(figsize=(8, 4), dpi=500)
        x = np.linspace(0, float(len(p.I)) / float(p.resolution), num=len(p.I))
        xcm = x / 0.000124
        Iabs = p.I.__abs__() / sum(p.I.__abs__()) * p.resolution
        if p.process == "absorption":
            plt.plot(xcm, Iabs, linewidth=1.5, label="PLE cross section", color="darkolivegreen")
            plt.title(name + "  PLE cross section")
            plt.xlim(plrange[0] / 0.000124, plrange[1] / 0.000124)
            plt.ylabel('Cross section($\hbar\omega$)')
        elif p.process == "emission":
            plt.plot(xcm, Iabs, linewidth=1.2, label="normalized PL lineshape", color="darkolivegreen")
            plt.title(name + "  PL lineshape ")
            plt.xlim(plrange[0] / 0.000124, plrange[1] / 0.000124)
            plt.ylabel(' $L_{em}(\hbar\omega)$')

        plt.xlabel('energy ($cm^{-1}$)')
        plt.ylim(0, max(Iabs) * 1.1)
        plt.vlines(p.EZPL / 0.000124, ymin=0, ymax=max(Iabs) * 1.1, linewidth=1, label='$E_{zpl}$=' + str(round(p.EZPL / 0.000124)) + "$cm^{-1}$", color='grey', linestyle="--")
        plt.legend()
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        plt.savefig("PLcm.png", dpi=500)

    if "PLnm" in image:
        plt.figure(figsize=(8, 4), dpi=500)
        wave = [1240.0 / q for q in x[int(0.5 * p.resolution):int(len(x) / 1.5)]]
        xnew = x[int(0.5 * p.resolution):int(len(x) / 1.5)]
        Iabs_new = Iabs[int(0.5 * p.resolution):int(len(x) / 1.5)]
        Iabswave = np.array([Iabs_new[i] * xnew[i]**2 for i in range(len(Iabs_new))])
        intergal = 0.0
        for i in range(len(Iabswave) - 1):
            intergal += Iabswave[i] * (wave[i + 1] - wave[i])
        print("intergal", intergal)
        Iabswave = Iabswave / abs(intergal)
        if p.process == "absorption":
            plt.plot(wave, Iabswave, label="PLE cross section", color="darkolivegreen")
            plt.title(name + "  PLE cross section")
            plt.xlim(1240.0 / plrange[1], 1240.0 / plrange[0])
            plt.ylabel('Cross section(nm)')
        elif p.process == "emission":
            plt.plot(wave, Iabswave, label=" PL intensity", color="darkolivegreen")
            plt.title(name + "   PL intensity")
            plt.xlim(1240.0 / plrange[1], 1240.0 / plrange[0])
            plt.ylabel(' PL intensity')
        data_ple = np.column_stack((wave, Iabswave, Iabswave / max(Iabswave)))
        np.savetxt('lineshape_nm.csv', data_ple, delimiter=',', header='wave lenth(nm),Iabs(a.u.),Iabs_norm(max=1)', comments='')
        plt.xlabel('wavelenth (nm) ')
        plt.ylim(0, max(Iabswave) * 1.1)
        plt.vlines(1240.0 / p.EZPL, ymin=0, ymax=max(Iabswave) * 1.1, label='$E_{zpl}$=' + str(round(1240.0 / p.EZPL, 1)) + "nm", color='grey', linestyle="--")
        plt.legend()
        plt.savefig("PLnm.png", dpi=500)