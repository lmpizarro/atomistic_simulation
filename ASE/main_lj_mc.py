#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division 

import os
import sys

pwd = os.getcwd()
mylibs = pwd + "/libs"
sys.path.append(mylibs)

import numpy as np

from ase.lattice.cubic import FaceCenteredCubic
import ase.data
from ase import Atoms

import ase.io

import lj_mc as ljmc
import lj_mix as ljmix

# init lennard jones interactions
elements = ["Ar", "Ne", "Kr", "Xe"]
interacciones = ljmix.Interacciones_LJ(1.2, 4.0, 2.3, elements[0])
interacciones.add_component(1.1, 4.5, 2.1, elements[1])
interacciones.add_component(1.0, 4.3, 2.0, elements[2])
interacciones.add_component(1.05, 4.4, 2.7, elements[3])

# init atoms systems
Nx=Ny=Nz=4
lc = 1.0 
atoms = FaceCenteredCubic(size=(Nx, Ny, Nz), symbol="Ne", pbc=True, latticeconstant=lc)
cs = [elements[0], elements[1],  elements[2], elements[3] ] * (Nx ** 3)
atoms.set_chemical_symbols(cs)		
atoms.info["elements"] = elements
atoms.info["pos_changed"] = 0
# fin de la inicializacion

# Associate the EMT potential with the atoms
atoms.set_calculator(ljmc.LennardJones_mc(interacciones))
print "energia inicial", atoms.get_potential_energy()/len(atoms)
cell = atoms.get_cell()
# ajuste para el minimo de energía potencial
f = 1.5 
atoms.positions *= f
atoms.set_cell(cell * f)
print "energia inicial", f, atoms.get_potential_energy()/len(atoms)



#atoms[0].position += 0.00001 
cs = atoms.get_chemical_symbols()
cs[0]="Ne"
atoms.set_chemical_symbols(cs)		
print "cambio de energia por cambio de elemento", atoms.get_potential_energy()/len(atoms)
atoms[0].position -= np.random.rand(3) 
atoms.info["pos_changed"] = 0
print "cambio de energia por cambio de posicion", atoms.get_potential_energy()/len(atoms)



ase.io.write('slab.xyz', atoms)
