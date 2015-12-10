#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import cubic


class Cubic(object):
    def __init__(self, lado_celda = 4, nx = 4, ny = 4, nz = 4, tipo="fcc"):
       self.Nperx = nx
       self.Npery = ny
       self.Nperz = nz

       self.lado_celda = lado_celda

       self.lado_cubo_x = self.lado_celda * self.Nperx
       self.lado_cubo_y = self.lado_celda * self.Npery
       self.lado_cubo_z = self.lado_celda * self.Nperz

       if (tipo == "fcc"):
           self.AtomsCell = 4
       elif (tipo =="simple"):
           self.AtomsCell = 1
       elif (tipo =="bcc"):
           self.AtomsCell = 2

       self.Natoms = self.AtomsCell * self.Nperx * self.Npery * self.Nperz
       self.Positions = np.zeros(3*self.Natoms, dtype='float64').reshape(3,self.Natoms)
       self.IndiceElementos = np.ones(self.Natoms, dtype='int')
       self.NombreElementos = ["C"] * self.Natoms

       if (tipo == "fcc"):
           self.Positions = cubic.posiciones_fcc(self.Positions, self.Nperx, 
                   self.Npery, self.Nperz, 3, self.Natoms) 
       elif (tipo =="simple"):
           self.Positions = cubic.posiciones_simple(self.Positions, self.Nperx, 
                   self.Npery, self.Nperz, 3, self.Natoms) 
       elif (tipo =="bcc"):
           self.Positions = cubic.posiciones_bcc(self.Positions, self.Nperx, 
                   self.Npery, self.Nperz, 3, self.Natoms) 

    def list_positions (self):
        print self.IndiceElementos
        print self.NombreElementos
        for i in range(np.size(self.Positions[0])):
            print  self.IndiceElementos[i], self.Positions[0][i], \
                  self.Positions[1][i], self.Positions[2][i]

    # representacion xyz de la celda básica
    def write_xyz(self, name):
        fo = open(name, 'w')
        fo.write(str(self.Natoms) + "\n")
        fo.write ("# xyz file\n")
        for i in range(np.size(self.Positions[0])):
            tmp = ("%s %f %f %f \n")%(self.IndiceElementos[i], self.Positions[0][i], \
                  self.Positions[1][i], self.Positions[2][i])
            fo.write( tmp )
        fo.close()

    # representacion vtf de la celda básica
    def write_vtf(self, name):
        fo = open(name, 'w')
        fo.write ("## vtf file\n")

        for i in range (self.Natoms):
            atom = ( " atom %d radius 1.0 type 1\n" )%(i)
            fo.write(atom)

        pbc = ("%s %f %f %f \n\n")%( "pbc", self.Nperx, self.Npery, self.Nperz)
        fo.write(pbc)
        fo.write (" timestep\n")
        for i in range(np.size(self.Positions[0])):
            tmp = ("   %f %f %f \n")%(self.Positions[0][i], \
                  self.Positions[1][i], self.Positions[2][i])
            fo.write( tmp )
        fo.close()

    def set_A4 (self):
        pass

    def set_A2B2(self):
        pass

    def set_A3B1(self):
        pass

def main ():
    fcc = Cubic(3, tipo="fcc")
    print fcc.Natoms
    print fcc.Nperx, fcc.Npery, fcc.Nperz
    print fcc.Positions
    print fcc.Positions.shape

    fcc.write_vtf("test.vtf")

    fcc.list_positions()


if __name__ == "__main__":
    main()
