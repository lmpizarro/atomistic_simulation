from __future__ import division 
import numpy as np

from ase.lattice.cubic import FaceCenteredCubic
import ase.data
from ase import Atoms

import ase.io

import lj_mc as ljmc
import lj_mix as ljmix

Nx=Ny=Nz=3
lc = 3/Nx
elements = ["Ar", "Ne", "Kr", "Xe"]

interacciones = ljmix.Interacciones_LJ(1.2, 4.0, 2.3, elements[0])
interacciones.add_component(1.1, 4.5, 2.1, elements[1])
interacciones.add_component(1.0, 4.3, 2.0, elements[2])
interacciones.add_component(1.05, 4.4, 2.7, elements[3])


atoms = FaceCenteredCubic(size=(Nx, Ny, Nz), symbol="Ne", pbc=True, latticeconstant=lc)
cs = [elements[0], elements[1],  elements[2], elements[3] ] * (Nx ** 3)
atoms.set_chemical_symbols(cs)		

# Associate the EMT potential with the atoms
atoms.set_calculator(ljmc.LennardJones_mc(interacciones))
print atoms.get_calculator().results['energy']

print atoms.get_potential_energy()/len(atoms)
#atoms[0].position = [1,1,1]
print atoms.get_potential_energy()/len(atoms)


#for i,a in enumerate(atoms.positions):
#  print i, a, cs[i]


ase.io.write('slab.xyz', atoms)

#print indices, offsets

print

#print nl.get_neighbors(3)
