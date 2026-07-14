import os
import numpy as np


class ConfigurationCoordinate:
    def read_poscar(self, i_path, l_get_sorted_symbols=False):
        from pymatgen.io.vasp.inputs import Poscar
        poscar = Poscar.from_file("{}".format(i_path))
        struct = poscar.structure
        if l_get_sorted_symbols:
            return struct, poscar.site_symbols
        else:
            return struct

    def Delta_Q(self, i_file, f_file, disp_range=None):
        struct_i, sorted_symbols = self.read_poscar(i_file, True)
        struct_f, sorted_symbols = self.read_poscar(f_file, True)
        delta_R = struct_f.frac_coords - struct_i.frac_coords
        delta_R = (delta_R + 0.5) % 1 - 0.5

        lattice = struct_i.lattice.matrix
        delta_R = np.dot(delta_R, lattice)

        masses = np.array([spc.atomic_mass for spc in struct_i.species])
        delta_Q2 = masses[:, None] * delta_R ** 2
        return np.sqrt(delta_Q2.sum())

    def get_init_fin(self, i_file, f_file, disp_range=None, output_dir='disp_dir'):
        from pymatgen.core.structure import Structure
        from pymatgen.io.vasp.inputs import Poscar
        if disp_range is None:
            disp_range = np.linspace(-1, 1, 11)
        struct_i, sorted_symbols = self.read_poscar(i_file, True)
        struct_f, sorted_symbols = self.read_poscar(f_file, True)
        delta_R = struct_f.frac_coords - struct_i.frac_coords
        delta_R = (delta_R + 0.5) % 1 - 0.5

        lattice = struct_i.lattice.matrix
        delta_R = np.dot(delta_R, lattice)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        masses = np.array([spc.atomic_mass for spc in struct_i.species])
        delta_Q2 = masses[:, None] * delta_R ** 2

        print('Delta_Q^2', np.sqrt(delta_Q2.sum()))

        for frac in disp_range:
            disp = frac * delta_R
            print(disp[0][0])
            struct = Structure(struct_i.lattice, struct_i.species,
                               struct_i.cart_coords + disp,
                               coords_are_cartesian=True)
            Poscar(struct).write_file('{0}/POSCAR_{1:03d}'.format(output_dir, int(np.rint(frac * 10))))