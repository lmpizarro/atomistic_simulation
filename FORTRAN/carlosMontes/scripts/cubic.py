#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import cubic

# http://www.petveturas.com/prog/scripts/create_crystal.py


class Cubic(object):

    def __init__(self, parametro_red=4, nx=4, ny=4, nz=4, tipo="fcc"):
        self.Nperx = nx
        self.Npery = ny
        self.Nperz = nz

        self.parametro_red = parametro_red

        self.lado_cubo_x = self.parametro_red * self.Nperx
        self.lado_cubo_y = self.parametro_red * self.Npery
        self.lado_cubo_z = self.parametro_red * self.Nperz

        if (tipo == "fcc"):
            self.AtomsCell = 4
        elif (tipo == "simple"):
            self.AtomsCell = 1
        elif (tipo == "bcc"):
            self.AtomsCell = 2

        self.Natoms = self.AtomsCell * self.Nperx * self.Npery * self.Nperz
        self.Positions = np.zeros(
            3 * self.Natoms,
            dtype='float64').reshape(
            3,
            self.Natoms)
        self.IndiceElementos = np.ones(self.Natoms, dtype='int')
        self.NombreElementos = ["C"] * self.Natoms

        if (tipo == "fcc"):
            self.Positions = self.generate_simple()
            f1 = self.generate_simple() + \
                np.array([0.5, 0.5, 0.0]).reshape(3, 1)
            f2 = self.generate_simple() + \
                np.array([0.5, 0.0, 0.5]).reshape(3, 1)
            f3 = self.generate_simple() + \
                np.array([0.0, 0.5, 0.5]).reshape(3, 1)
            self.Positions = np.hstack((self.Positions, f1))
            self.Positions = np.hstack((self.Positions, f2))
            self.Positions = np.hstack((self.Positions, f3))
        elif (tipo == "simple"):
            self.Positions = self.generate_simple()
        elif (tipo == "bcc"):
            self.Positions = self.generate_simple()
            f1 = self.generate_simple() + \
                np.array([0.5, 0.5, 0.5]).reshape(3, 1)

        self.Positions = self.parametro_red * self.Positions

    def list_positions(self):
        print self.IndiceElementos
        print self.NombreElementos
        for i in range(np.size(self.Positions[0])):
            print  self.IndiceElementos[i], self.Positions[0][i], \
                self.Positions[1][i], self.Positions[2][i]

    def generate_simple(self):
        p = np.zeros(
            3 *
            self.Nperx *
            self.Npery *
            self.Nperz).reshape(
            3,
            self.Nperx *
            self.Npery *
            self.Nperz)
        igr = 0
        for i in range(self.Nperx):
            for j in range(self.Npery):
                for k in range(self.Nperz):
                    p[0][igr] = i
                    p[1][igr] = j
                    p[2][igr] = k
                    igr += 1
        return p

        # print p.shape

    # representacion xyz de la celda básica
    def write_xyz(self, name):
        fo = open(name, 'w')
        fo.write(str(self.Natoms) + "\n")
        fo.write("# xyz file\n")
        for i in range(np.size(self.Positions[0])):
            tmp = ("%s %f %f %f \n") % (self.IndiceElementos[i], self.Positions[
                0][i], self.Positions[1][i], self.Positions[2][i])
            fo.write(tmp)
        fo.close()

    # representacion vtf de la celda básica
    def write_vtf(self, name):
        fo = open(name, 'w')
        fo.write("## vtf file\n")

        for i in range(self.Natoms):
            atom = (" atom %d radius %f type %d\n") % (i,
                    self.radios[self.IndiceElementos[i] - 1],\
                    self.IndiceElementos[i])
            fo.write(atom)

        pbc = ("%s %f %f %f \n\n") % (
            "pbc", self.Nperx, self.Npery, self.Nperz)
        fo.write(pbc)
        fo.write(" timestep\n")
        for i in range(np.size(self.Positions[0])):
            tmp = ("%f %f %f \n") % (
                self.Positions[0][i], self.Positions[1][i], self.Positions[2][i])
            fo.write(tmp)
        fo.close()

    # escribe la configuracion de posiciones
    # atómicas
    def write_config(self, name):
        format1 = "%14.8f %14.8f %14.8f %10d\n"
        format2 = "%10d %10d %10d %10d %14.8f\n"
        fo = open(name, 'a')
        fo.write((format2) % (self.Natoms, self.Nperx, self.Npery, self.Nperz,\
            self.parametro_red))
        for i in range(np.size(self.Positions[0])):
            tmp = (format1) % (self.Positions[0][i],
                              self.Positions[1][i], self.Positions[2][i],
                              self.IndiceElementos[i])
            fo.write(tmp)
        fo.close()

    def set_A4(self):
        pass

    def set_A2B2(self):
        self.radios=np.array([.8, 1.0])
        for i in range(self.Natoms / 2):
            self.IndiceElementos[i] = 2

    def set_A3B1(self):
        pass


def main():
    fcc = Cubic(2, nx=4, ny=4, nz=4, tipo="fcc")
    print fcc.Natoms
    print fcc.Nperx, fcc.Npery, fcc.Nperz
    print fcc.Positions
    print fcc.Positions.shape

    fcc.set_A2B2()

    fcc.write_xyz("test.xyz")
    fcc.write_vtf("test.vtf")
    fcc.write_config("config.par")

    fcc.list_positions()


if __name__ == "__main__":
    main()
