#!/usr/bin/python
# -*- coding: utf-8 -*-
from types import *

import numpy as np
import cubic
import sys

import ase.data as aseData
from ase import Atoms
from ase import Atom

# http://www.petveturas.com/prog/scripts/create_crystal.py

symbol_list = aseData.chemical_symbols

class Cubic(object):

    def __init__(self, parametro_red=4, nx=4, ny=4, nz=4, tipo="fcc"):

        try: 
            assert type(parametro_red) is FloatType, "parametro_red must be float"
            assert (parametro_red > 0), "Must be greater than zero!!"

            assert type(nx) is IntType, "nx must be integer"
            assert (nx >= 1), "Must be greater or equal one (1) !!"

            assert type(ny) is IntType, "ny must be integer"
            assert (ny >= 1), "Must be greater or equal one (1) !!"

            assert type(nz) is IntType, "nz must be integer"
            assert (nz >= 1), "Must be greater or equal one (1) !!"

            assert type(tipo) is StringType, "tipo must be string"
            assert (tipo == "fcc" or tipo == "simple" or tipo == "bcc") 

        except AssertionError as e:
            print ("%s %s ") %("Sistema mal configurado", e)
            sys.exit(2)


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
        self.Positions = np.zeros(3 * self.Natoms, dtype='float64').reshape(
            3, self.Natoms)
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

    def set_parametro_red(self, parametro_red):

        try:
            assert type(parametro_red) is FloatType, "parametro_red must be float"
            assert (parametro_red > 0), "Must be greater than zero!!"
        except AssertionError as e:
            print ("%s %s ") %("Bad net parameter: ", e)
            sys.exit(2)

        tmp_pr = self.parametro_red

        self.parametro_red = parametro_red

        self.lado_cubo_x = self.parametro_red * self.Nperx
        self.lado_cubo_y = self.parametro_red * self.Npery
        self.lado_cubo_z = self.parametro_red * self.Nperz

        self.Positions = self.parametro_red * self.Positions/tmp_pr 


    def list_positions(self):
        # print self.IndiceElementos
        # print self.NombreElementos
        for i in range(np.size(self.Positions[0])):
            print  self.IndiceElementos[i], self.Positions[0][i], \
                self.Positions[1][i], self.Positions[2][i]

    def generate_simple(self):
        p = np.zeros(3 * self.Nperx * self.Npery * self.Nperz).reshape(
            3, self.Nperx * self.Npery * self.Nperz)
        igr = 0
        for i in range(self.Nperx):
            for j in range(self.Npery):
                for k in range(self.Nperz):
                    p[0][igr] = i
                    p[1][igr] = j
                    p[2][igr] = k
                    igr += 1
        return p


    # representacion xyz de la celda básica
    def write_xyz(self, name):

        try:
            assert type(name) is StringType, "name must be string"
            end_name = name[-4:]
            assert end_name == ".xyz", "config file must end with .par"
        except AssertionError as e:
            print ("%s %s ") %(" Bad file name: ", e)
            sys.exit(2)

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
        try:
            assert type(name) is StringType, "name must be string"
            end_name = name[-4:]
            assert end_name == ".vtf", "config file must end with .par"
        except AssertionError as e:
            print ("%s %s ") %(" Bad file name: ", e)
            sys.exit(2)


        fo = open(name, 'w')
        fo.write("## vtf file\n")

        for i in range(self.Natoms):
            atom = (" atom %d radius %f type %d\n") % (i,
                                                       self.covalent_radii[
                                                           self.IndiceElementos[i] - 1],
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

        try:
            assert type(name) is StringType, "name must be string"
            end_name = name[-4:]
            assert end_name == ".par", "config file must end with .par"
        except AssertionError as e:
            print ("%s %s ") %(" Bad file name: ", e)
            sys.exit(2)

        format1 = "%14.8f %14.8f %14.8f %10d\n"
        format2 = "%10d %10d %10d %10d %14.8f\n"
        fo = open(name, 'a')
        fo.write((format2) % (self.Natoms, self.Nperx, self.Npery, self.Nperz,
                              self.parametro_red))
        for i in range(np.size(self.Positions[0])):
            tmp = (format1) % (self.Positions[0][i],
                               self.Positions[1][i], self.Positions[2][i],
                               self.IndiceElementos[i])
            fo.write(tmp)
        fo.close()

    def cpc(self):
        #tmp = gLado_caja*floor(gR/gLado_caja)
        #gR = gR - tmp
        pass
        tmpx = self.lado_cubo_x * \
            np.floor(self.Positions[0, :] / self.lado_cubo_x)
        tmpy = self.lado_cubo_x * \
            np.floor(self.Positions[1, :] / self.lado_cubo_y)
        tmpz = self.lado_cubo_x * \
            np.floor(self.Positions[2, :] / self.lado_cubo_z)

        tmp = np.zeros(3 * self.Natoms).reshape(3, self.Natoms)
        tmp[0, :] = tmpx
        tmp[1, :] = tmpy
        tmp[2, :] = tmpz

        self.Positions = self.Positions - tmp

    def noisify(self, amplitude):

        try:
            assert (amplitude > 0), "Amplitude must be greater than 0"
            assert type(amplitude) is  FloatType, "Amplitude must be float"
        except AssertionError as e:
            print ("%s %s ") %("Wrong Amplitude Noise: ", e)
            sys.exit(2)

        tmp = amplitude * \
            np.random.rand(self.Natoms * 3). reshape(3, self.Natoms)
        self.Positions = self.Positions + tmp
        self.cpc()

    def set_A4(self):
        pass

    # por ahora valido para fcc
    def set_A2B2(self, a1, a2):
        self.radios = np.array([.8, 1.0])
        for i in range(self.Natoms / 2):
            self.IndiceElementos[i] = 2

    
    class AtomNameError(Exception):
            def __init__(self, value):
                self.value = value
            def __str__(self):
                return repr(self.value)

    def set_ase_atoms_rep(self):
        self.ase_atoms = Atoms(pbc=True)
        for i in range(np.size(self.Positions[0])):
            position = [self.Positions[0,i], self.Positions[1,i], self.Positions[2,i]] 
            symb = self.symbols[self.IndiceElementos[i] - 1]
            a = Atom(symb)
            a.position = position 
            self.ase_atoms.append(a)


    # por ahora valido para fcc
    def set_A1B1C1D1(self, a1, a2, a3, a4):

        try:
            assert type(a1) is  StringType, "atom 1  must be string"
            assert type(a2) is  StringType, "atom 2  must be string"
            assert type(a3) is  StringType, "atom 3  must be string"
            assert type(a4) is  StringType, "atom 4  must be string"
        except AssertionError as e:
            print ("%s %s ") %("Wrong atoms symbols", e)
            raise self.AtomNameError(e)

        try:
            assert a1 in symbol_list, "atom 1  not in symbol list"
            assert a2 in symbol_list, "atom 2  not in symbol list"
            assert a3 in symbol_list, "atom 3  not in symbol list"
            assert a4 in symbol_list, "atom 4  not in symbol list"
        except AssertionError as e:
            print ("%s %s ") %("Wrong atoms symbols", e)
            raise self.AtomNameError(e)


        self.symbols = [a1, a2, a3, a4]
        self.atomic_numbers = []
        self.atomic_masses  = []
        self.covalent_radii = []

        [self.atomic_numbers.append(aseData.atomic_numbers[s]) for s in self.symbols]
        [self.atomic_masses.append(aseData.atomic_masses[n]) for n in self.atomic_numbers]
        [self.covalent_radii.append(aseData.covalent_radii[n]) for n in self.atomic_numbers]

        for i in range(self.Natoms / 4):
            self.IndiceElementos[i] = 2
        for i in range(self.Natoms/4,self.Natoms / 2):
            self.IndiceElementos[i] = 3
        for i in range(self.Natoms/2,3*self.Natoms / 4):
            self.IndiceElementos[i] = 4

        self.set_ase_atoms_rep()    

    def set_A3B1(self):
        pass


def main():
    fcc = Cubic(2.0, nx=4, ny=4, nz=4, tipo="fcc")
    #print fcc.Natoms
    #print fcc.Nperx, fcc.Npery, fcc.Nperz
    #print fcc.Positions
    #print fcc.Positions.shape

    fcc.set_A1B1C1D1("Ar","Ne","Kr","Xe")
    #fcc.noisify(2.0)

    fcc.write_xyz("test.xyz")
    fcc.write_vtf("test.vtf")
    fcc.write_config("config.par")

    #fcc.list_positions()
    print fcc.ase_atoms.positions

if __name__ == "__main__":
    main()
