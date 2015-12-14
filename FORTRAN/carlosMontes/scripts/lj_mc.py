#!/usr/bin/python
# -*- coding: utf-8 -*-


from numba import jit
import numpy as np


class LJ(object):

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
                #print (rij)
                if ( d2 < Rc2):
                    pass
                    id2 = 1. / d2            #inverse squared distance
                    id6 = id2 * id2 * id2    #inverse sixth distance
                    id12 = id6 * id6         #inverse twelvth distance
                    PEnergy = PEnergy + 4. * (id12 - id6) + Shift
                pass
        return(PEnergy)

    def delta_energy(self)
        DEnergy = 0.0
        return(DEnergy)
