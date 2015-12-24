from __future__ import division
import numpy as np
from ase.calculators.neighborlist import NeighborList 
from ase.calculators.calculator import Calculator, all_changes 
import types
#import lj_mix
import sys

class LJ_MIX_Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class LennardJones_mc(Calculator): 
    implemented_properties = ['energy']

    default_parameters = {'epsilon': 1.0,
                          'sigma': 1.0,
                          'rc': None}

    nolabel = True

    def __init__(self, lj_mix, **kwargs):
        Calculator.__init__(self, **kwargs)
	energy = 10

        self.results['energy'] = energy

	a = ("%s")%(type(lj_mix)) 
	try:
	    assert a == "<class 'lj_mix.Interacciones_LJ'>"
        except AssertionError as e:
	    print ("%s %s ") %("It isn't a lj_mix instance", e)
            raise self.LJ_MIX_Error(e)

        self.lj_mix = lj_mix
        self.lj_mix.calculate_matrices()
	print self.lj_mix


    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=['positions']):
	print "calculando"
        Calculator.calculate(self, atoms, properties, system_changes)


        natoms = len(self.atoms)

        sigma = self.parameters.sigma
        epsilon = self.parameters.epsilon
        rc = self.parameters.rc
        if rc is None:
            rc = 3 * sigma
        
        if 'numbers' in system_changes:
            self.nl = NeighborList([rc / 2] * natoms, self_interaction=False)

        self.nl.update(self.atoms)
        
        positions = self.atoms.positions
        cell = self.atoms.cell
        
        e0 = 4 * epsilon * ((sigma / rc)**12 - (sigma / rc)**6)
        
        energy = 0.0
        
        
        self.results['energy'] = energy


