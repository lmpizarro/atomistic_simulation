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
#σ = 3.4 × 10-10 m and ε = 1.65 × 10-21 joule
# 
elements = ["Ar", "Ne", "Kr", "Xe"]
interacciones = ljmix.Interacciones_LJ(0.01, 3.4, 2.3, elements[0])
interacciones.add_component(0.01, 4.5, 2.1, elements[1])
interacciones.add_component(0.03, 4.3, 2.0, elements[2])
interacciones.add_component(0.02, 4.4, 2.7, elements[3])
interacciones.calculate_matrices()

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
calculator = ljmc.LennardJones_mc(interacciones)
atoms.set_calculator(calculator)
print "energia inicial", atoms.get_potential_energy()/len(atoms)
cell = atoms.get_cell()
# ajuste para el minimo de energía potencial
f = 1.41
atoms.positions *= f
atoms.set_cell(cell * f)
Einic = atoms.get_potential_energy()/len(atoms)
print "energia inicial", f, Einic



#atoms[0].position += 0.00001 
cs = atoms.get_chemical_symbols()
cs[0]="Ar"
atoms.set_chemical_symbols(cs)		
print "cambio de energia por cambio de elemento", atoms.get_potential_energy()/len(atoms)

aceptados = 0
Energy = 0
kT = 0.005103
beta = 1.0 /kT
Energy = Einic
for i in range (100000):
  indx_change = np.random.randint(0,Nx*Nz*Ny*4)
 
  atoms.info["pos_changed"] = indx_change
  atoms.info["pos_changed_value"] = np.copy(atoms[indx_change].position)
  p = atoms[indx_change].position
  p_cpc = calculator.energy_lj.cpc(p)
  atoms[indx_change].position = p_cpc
 
  einic =  atoms.get_potential_energy()/len(atoms)

  
  atoms[indx_change].position += 0.02 * lc * f * (np.random.rand(3) - 0.5)

  p = atoms[indx_change].position
  p_cpc = calculator.energy_lj.cpc(p)
  atoms.positions[indx_change] = p_cpc
  eifin = atoms.get_potential_energy()/len(atoms)


  deltaE = eifin - einic
  	

  if deltaE < 0:
     aceptados += 1
     Energy += deltaE 
     print aceptados, i, deltaE, Energy
  else:
     if np.random.rand() < np.exp(-deltaE*beta):
        aceptados += 1
        Energy += deltaE 
        print aceptados, i, deltaE, Energy


#print atoms.positions

ase.io.write('slab.xyz', atoms)
