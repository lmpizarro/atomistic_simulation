#!/usr/bin/python
# -*- coding: utf-8 -*-
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


	a = ("%s")%(type(lj_mix)) 
	try:
	    assert a == "<class 'lj_mix.Interacciones_LJ'>"
        except AssertionError as e:
	    print ("%s %s ") %("It isn't a lj_mix instance", e)
            raise self.LJ_MIX_Error(e)

        self.lj_mix = lj_mix
        self.lj_mix.calculate_matrices()
	print self.lj_mix

        self.results['energy'] = energy

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

	print system_changes    

        self.nl.update(self.atoms)
        
        positions = self.atoms.positions
        cell = self.atoms.cell
        
        e0 = 4 * epsilon * ((sigma / rc)**12 - (sigma / rc)**6)
        
        energy = 0.0
        
        
        self.results['energy'] = energy

from numba import jit
import numpy as np


class EnergyLJ(object):

    def __init__(self, components, system):
        self.componentes = components
        self.system = system
        pass

    @jit
    def cpc(self, r):
        tmpx = self.system.lado_cubo_x * \
            np.floor(r[0] / self.system.lado_cubo_x)
        tmpy = self.system.lado_cubo_x * \
            np.floor(r[1] / self.system.lado_cubo_y)
        tmpz = self.system.lado_cubo_x * \
            np.floor(r[2] / self.system.lado_cubo_z)
        
        tmp = np.array([tmpx, tmpy, tmpz])

        r = r - tmp

        return (r)

    @jit
    def energy(self):
        PEnergy = 0.0
        for i in range(self.system.Natoms - 1):
            Posi = self.system.Positions[:, i]
            for j in range(i + 1, self.system.Natoms):
                rij = self.system.Positions[:,j] - Posi
                Shift = self.componentes.potCut[self.system.IndiceElementos[i],
                                                self.system.IndiceElementos[j]]
                Epsilon = self.componentes.epsilon[
                    self.system.IndiceElementos[i],
                    self.system.IndiceElementos[j]]
                Sigma = self.componentes.sigma[self.system.IndiceElementos[i],
                                               self.system.IndiceElementos[j]]
                Rc2 = self.componentes.rc2[self.system.IndiceElementos[i],
                                           self.system.IndiceElementos[j]]
                rij = self.cpc(rij)
                d2 = Sigma * np.dot(rij, rij)
                if ( d2 < Rc2):
                    pass
                    id2 = 1. / d2            #inverse squared distance
                    id6 = id2 * id2 * id2    #inverse sixth distance
                    id12 = id6 * id6         #inverse twelvth distance
                    PEnergy = PEnergy + (4. * (id12 - id6))*Epsilon + Shift
                pass
        return(PEnergy)

    @jit
    def kernel_e (self, i, j):
        Posi = self.system.Positions[:, i]
        rij = self.system.Positions[:,j] - Posi
        Shift = self.componentes.potCut[self.system.IndiceElementos[i],
                                                self.system.IndiceElementos[j]]
        Epsilon = self.componentes.epsilon[ self.system.IndiceElementos[i],
                    self.system.IndiceElementos[j]]
        Sigma = self.componentes.sigma[self.system.IndiceElementos[i],
                                               self.system.IndiceElementos[j]]
        Rc2 = self.componentes.rc2[self.system.IndiceElementos[i],
                                           self.system.IndiceElementos[j]]
        rij = self.cpc(rij)
        d2 = Sigma * np.dot(rij, rij)
        if ( d2 < Rc2):
            id2 = 1. / d2            #inverse squared distance
            id6 = id2 * id2 * id2    #inverse sixth distance
            id12 = id6 * id6         #inverse twelvth distance
            return ((4. * (id12 - id6))*Epsilon + Shift)
        else:
            return 0


    @jit
    def energy_pos (self, ii):
        e_pos = 0.0
        for j in range(ii):
            e_pos = e_pos + self.kernel_e(ii, j)

        for j in range(ii + 1, self.system.Natoms):
            e_pos = e_pos + self.kernel_e(ii, j)

        return e_pos


    @jit
    def energy_IJ(self, ii, jj):
        PEnergy = 0.0

        # hacer un calculo para ii
        # átomo en la posicion ii sin cambiar

        E1 = self.energy_pos(ii)

        # hacer un calculo para ii
        # cambiando la identidad por el de jj
        idxi = self.system.IndiceElementos[ii] 
        idxj = self.system.IndiceElementos[jj] 
        self.system.IndiceElementos[ii] = idxj 
        self.system.IndiceElementos[jj] = idxi
        E2 = self.energy_pos(ii)
        
        # vuelve los identificadores a su estado
        # original
        #self.system.IndiceElementos[ii] = idxi 
        #self.system.IndiceElementos[jj] = idxj

        return (E2 - E1)

    @jit
    def distance_true(self, i, j):
        Sigma = self.componentes.sigma[self.system.IndiceElementos[i],
                                               self.system.IndiceElementos[j]]
        Rc2 = self.componentes.rc2[self.system.IndiceElementos[i],
                                           self.system.IndiceElementos[j]]

        Posi = self.system.Positions[:, i]
        rij = self.system.Positions[:,j] - Posi
        rij = self.cpc(rij)
        d2 = Sigma * np.dot(rij, rij)

        return d2 < Rc2


