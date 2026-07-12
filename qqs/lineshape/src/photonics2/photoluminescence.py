import numpy as np
from numpy import fft
from photonics2.xyz import XYZ
from photonics2.configuration_coordinate import ConfigurationCoordinate
from scipy import constants as constant


def realfft(func, emax):
    return fft.fft(func) / float(len(func)) * emax


def realifft(func, emax):
    return fft.ifft(func) * emax


def gaussian(omega, omega_k, sigma):
    return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(omega - omega_k)**2 / sigma**2 / 2)


def S_omega_vectorized(omega_set, sigma, frequencies, S, skipmodes, jtmodes):
    mask = np.ones(len(frequencies), dtype=bool)
    mask[:skipmodes + 1] = False
    if jtmodes:
        mask[list(jtmodes)] = False

    freqs = frequencies[mask]
    Sk = S[mask]
    sigma_k = np.where(freqs < 0.005, sigma, 0.3 * sigma)

    diff = omega_set[:, np.newaxis] - freqs[np.newaxis, :]
    weight = Sk[np.newaxis, :] / (np.sqrt(2 * np.pi) * sigma_k[np.newaxis, :])
    result = weight * np.exp(-diff**2 / (2 * sigma_k[np.newaxis, :]**2))
    return np.sum(result, axis=1)


def C_omega_vectorized(omega_set, sigma, frequencies, S, skipmodes, temperature, jtmodes):
    mask = np.ones(len(frequencies), dtype=bool)
    mask[:skipmodes + 1] = False
    if jtmodes:
        mask[list(jtmodes)] = False

    freqs = frequencies[mask]
    Sk = S[mask]

    hbar = constant.physical_constants["Planck constant over 2 pi"][0]
    kB = constant.physical_constants["Boltzmann constant"][0]
    hbar_eV = constant.physical_constants["Planck constant over 2 pi in eV s"][0]
    a = hbar * freqs / temperature / kB / hbar_eV
    occupation = 1 / (np.exp(a) - 1)

    diff = omega_set[:, np.newaxis] - freqs[np.newaxis, :]
    weight = occupation[np.newaxis, :] * Sk[np.newaxis, :] / (np.sqrt(2 * np.pi) * sigma)
    result = weight * np.exp(-diff**2 / (2 * sigma**2))
    return np.sum(result, axis=1)


class Photoluminescence:
    _UNIT_THZ_TO_EV = 0.004135665538536
    _UNIT_SIESTA = 0.727445665

    def phonopy_read_modes(self):
        modes = np.zeros((self.numModes, self.numAtoms, 3))
        try:
            band = open(self.path, 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            raise

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            band.readline()
            band.readline()
            for a in range(self.numAtoms):
                band.readline()
                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][0] = float(parts[2])
                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][1] = float(parts[2])
                line = band.readline().replace(",", "")
                parts = line.strip().split()
                modes[i][a][2] = float(parts[2])

        band.close()
        return modes

    def phonopy_read_frequencies(self):
        frequencies = np.zeros(self.numModes)
        try:
            band = open(self.path, 'r')
        except OSError:
            print("Could not open/read file: band.yaml")
            raise

        for line in band:
            if "  band:" in line:
                break

        for i in range(self.numModes):
            band.readline()
            line = band.readline()
            parts = line.strip().split()
            frequencies[i] = float(parts[1])
            band.readline()
            for a in range(self.numAtoms):
                band.readline()
                band.readline()
                band.readline()
                band.readline()

        band.close()
        return frequencies

    def fold(self, d):
        nml = [round(a) for a in np.linalg.solve(self.lattice_vector.T, d)]
        return d - np.matmul(self.lattice_vector.T, nml)

    def get_phonon(self, path):
        with open(path, 'r') as f:
            natom = None
            masses = []

            for line in f:
                if line.startswith('natom:'):
                    natom = int(line.split(':')[1])
                elif line.strip() == 'points:':
                    for _ in range(natom):
                        next(f)
                        next(f)
                        line = next(f)
                        masses.append(float(line.split(':')[1]))
                elif "  band:" in line:
                    break

            numModes = 3 * natom
            m = np.array(masses) * constant.physical_constants["atomic mass constant"][0]
            frequencies = np.zeros(numModes)
            modes = np.zeros((numModes, natom, 3))

            for i in range(numModes):
                f.readline()
                line = f.readline()
                frequencies[i] = float(line.split(':')[1])
                f.readline()

                for a in range(natom):
                    f.readline()
                    line = f.readline().replace(',', '')
                    modes[i, a, 0] = float(line.strip().split()[2])
                    line = f.readline().replace(',', '')
                    modes[i, a, 1] = float(line.strip().split()[2])
                    line = f.readline().replace(',', '')
                    modes[i, a, 2] = float(line.strip().split()[2])

        sqrt_m = np.sqrt(m)
        modes = modes / sqrt_m[np.newaxis, :, np.newaxis]
        norms = np.linalg.norm(modes.reshape(numModes, -1), axis=1)
        norms[norms == 0] = 1
        modes = modes / norms[:, np.newaxis, np.newaxis]

        return m, modes, frequencies

    def read_grd_ex_pos(self, str_g, str_e, shiftVector, n_defect, defect_range, path):
        if '.xyz' in str_g:
            self.g = XYZ(str_g).coordinates
            self.e = XYZ(str_e).coordinates
        else:
            cc = ConfigurationCoordinate()
            self.g = cc.read_poscar(str_g)
            self.e = cc.read_poscar(str_e)
            self.lattice_vector = self.g.lattice.matrix

            self.g.translate_sites(range(len(self.g.frac_coords)), shiftVector, frac_coords=False)
            self.e.translate_sites(range(len(self.e.frac_coords)), shiftVector, frac_coords=False)

            lg = self.g.lattice
            le = self.e.lattice
            self.g = lg.get_cartesian_coords(self.g.frac_coords)
            self.e = le.get_cartesian_coords(self.e.frac_coords)

        self.numAtoms = len(self.g)
        D_R = []
        for i in range(self.numAtoms):
            D_R.append(self.fold(self.e[i] - self.g[i]).tolist())
        self.D_R = np.array(D_R)

        self.m, self.Modes, self.frequencies = self.get_phonon(path)

        if n_defect > 0:
            n_defect = n_defect - 1
            sumd = np.array([0.0, 0.0, 0.0])
            total_weight = 0.0
            for i in range(self.numAtoms):
                if np.dot(self.e[i] - self.e[n_defect], self.e[i] - self.e[n_defect]) > defect_range:
                    weight = self.m[i]
                    sumd += self.D_R[i] * weight
                    total_weight += weight
            sumd /= total_weight
            for i in range(self.numAtoms):
                self.D_R[i] -= sumd

    def el_ph(self, **parameter):
        self.temperature = parameter.get("temperature", 0)
        self.delta_width = parameter.get("delta_width", 6e-3)
        self.jtmodes = parameter.get("jtmodes", [])
        self.omega_set = np.linspace(
            self.min_energy, self.max_energy,
            int((self.max_energy - self.min_energy) * self.resolution))

        self.S_omega = S_omega_vectorized(
            self.omega_set, self.delta_width, self.frequencies,
            self.S, self.skipmodes, self.jtmodes)

        if self.temperature > 1e-2:
            self.C_omega = C_omega_vectorized(
                self.omega_set, self.delta_width, self.frequencies,
                self.S, self.skipmodes, self.temperature, self.jtmodes)
        else:
            self.C_omega = self.S_omega

        with open("C_omega.data", "w") as f:
            for a in self.C_omega:
                f.write(str(a) + "\n")

    def PL(self, **parameter):
        self.EZPL = parameter.get("EZPL", 0)
        self.gamma = parameter.get("gamma", 1)
        self.process = parameter.get("process", "emission")
        SHR = parameter.get("SHR", 0)

        r = 1 / self.resolution
        if "emission" in self.process:
            St = realfft(self.S_omega, self.max_energy - self.min_energy)
            Ct = realfft(self.C_omega, self.max_energy - self.min_energy)
            St = fft.fftshift(St)
            Ct = fft.fftshift(Ct)
        elif "absorption" in self.process:
            St = realifft(self.S_omega, self.max_energy - self.min_energy)
            Ct = realifft(self.C_omega, self.max_energy - self.min_energy)
            St = fft.fftshift(St)
            Ct = fft.fftshift(Ct)
        else:
            print("process input error, input:", self.process)
            raise ValueError(f"Unknown process: {self.process}")

        n = len(St)
        # Apply gamma broadening as a Lorentzian convolution in the
        # energy domain instead of exp(-gamma*|t|) in the time domain.
        # This decouples line-shape width from `resolution`. See
        # tools/RESOLUTION_NOTES.md for analysis.
        if self.temperature > 1e-2:
            G = np.exp(St + Ct + Ct[::-1] - St[0] - 2 * Ct[0] - SHR)
        else:
            G = np.exp(St - St[0] - SHR)
        A = fft.ifft(G) * n / (self.max_energy - self.min_energy)
        omega = np.arange(n) * r - (n // 2) * r
        kernel = (self.gamma / np.pi) / (omega**2 + self.gamma**2)
        kernel = kernel / kernel.sum()
        A = np.convolve(A, kernel, mode="same")
        shift = int((self.EZPL - self.min_energy) * self.resolution)
        idx = (shift - np.arange(n)) % n
        self.A = A[idx]

        return self.A

    def PLA(self):
        r = 1 / self.resolution
        t = r * (np.arange(len(self.A)) + self.min_energy * self.resolution)
        if "emission" in self.process:
            self.I = self.A * (t * r)**3
        elif "absorption" in self.process:
            self.I = self.A * (t * r)
        else:
            print("process input error, input:", self.process)
            raise ValueError(f"Unknown process: {self.process}")

        return self.I

    def __init__(self, path, method, **parameter):
        self.path = path
        self.method = method
        self.resolution = parameter.get("resolution", 1000)

        self.read_grd_ex_pos(
            parameter.get("POSCAR_GRD"),
            parameter.get("POSCAR_EX"),
            parameter.get('shift_vector=', [0.0, 0.0, 0.00]),
            parameter.get("n_defect", 0),
            parameter.get("defect_range", -1.0),
            path)

        if parameter.get('m', 0) != 0:
            self.numModes = parameter.get('nmodes', 3 * self.numAtoms)
            self.Modes = self.phonopy_read_modes()
            self.frequencies = self.phonopy_read_frequencies()
            self.m = parameter.get('m')
        else:
            self.m, self.Modes, self.frequencies = self.get_phonon(path)
            self.numModes = len(self.Modes)

        self.skipmodes = -1
        for i in range(self.numModes):
            if "vasp" in method:
                self.frequencies[i] = self.frequencies[i] / 1000
            elif "phonopy-siesta" in method:
                self.frequencies[i] = self.frequencies[i] * self._UNIT_THZ_TO_EV * self._UNIT_SIESTA
            elif "phonopy" in method:
                self.frequencies[i] = self.frequencies[i] * self._UNIT_THZ_TO_EV
            if self.frequencies[i] <= 0.00:
                self.skipmodes = i

        print(" modes energy < 0", self.skipmodes + 1)

    def HuangRhyes(self):
        self.HuangRhyes = 0
        self.Delta_R = 0
        self.Delta_Q = 0

        atomunit = np.sqrt(constant.physical_constants["atomic mass constant"][0]) / 1e10

        participation = np.sum(self.Modes**2, axis=2)
        self.IPR = 1.0 / np.sum(participation**2, axis=1)

        sqrt_m = np.sqrt(self.m)
        self.q = np.sum(
            sqrt_m[np.newaxis, :, np.newaxis] * self.D_R[np.newaxis, :, :] * self.Modes,
            axis=(1, 2)) * 1e-10

        self.q_amu = self.q / atomunit

        hbar = constant.physical_constants["Planck constant over 2 pi"][0]
        hbar_eV = constant.physical_constants["Planck constant over 2 pi in eV s"][0]
        self.S = self.frequencies * self.q**2 / 2 / (hbar * hbar_eV)
        # Zero out negative-frequency modes (pyphotonics clips freqs<0 to 0).
        self.S = np.where(self.frequencies < 0, 0.0, self.S)
        self.HuangRhyes = np.sum(self.S)

        with open("D(e-g).data", "w") as f:
            for a in range(self.numAtoms):
                dr_sq = np.sum(self.D_R[a]**2)
                self.Delta_R += dr_sq
                self.Delta_Q += dr_sq * self.m[a]
                f.write(str(self.D_R[a]) + "\t\n")

        self.Delta_R = self.Delta_R ** 0.5
        self.Delta_Q = (self.Delta_Q / constant.physical_constants["atomic mass constant"][0]) ** 0.5

        self.max_energy = 5.0
        self.min_energy = -0.0

        return self.HuangRhyes

    def print_table(self, k, visualmode):
        atomunit = np.sqrt(constant.physical_constants["atomic mass constant"][0]) / 1e10
        with open("main_modes.data", "w") as f, open("modes.data", "w") as f2:
            f.write("Modes NO.\t DeltaQ(A^2)\t IPR\t energy(meV)\t partial HuangRhyes\n ")
            f2.write("Modes NO.\t DeltaQ(A^2)\t IPR\t energy(meV)\t partial HuangRhyes\n ")
            for i in range(self.numModes):
                f2.write(f" {i + 1}\t{self.q[i] / atomunit}\t{self.IPR[i]} \t{self.frequencies[i] * 1000} \t{self.S[i]}\n")
                if self.S[i] > k and self.frequencies[i] > 0.04 or i in visualmode:
                    print(f"mode No.{i + 1} q={self.q[i] / atomunit} IPR={round(self.IPR[i], 2)} energy={round(self.frequencies[i] * 1000, 3)}meV  Sk={self.S[i]}")
                    f.write(f" {i + 1}\t{self.q[i] / atomunit}\t{self.IPR[i]} \t{self.frequencies[i]} \t{self.S[i]}\n")
