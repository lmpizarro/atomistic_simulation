#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys

pwd = os.getcwd()
mylibs = pwd + "/libs"
sys.path.append(mylibs)


import lj_mc as LJ
import paramCorrida as PC
import lj_mix as PM
import cubic as CB
import numpy as np


import time


def main():
    # inicializo los componentes del cálculo
    components = PM.Interacciones_LJ(1.1, 4.0, 2.3, "componente1")
    components.add_component(1.1, 4.5, 2.1, "componente2")
    components.add_component(1.0, 4.3, 2.0, "componente2")
    components.calculate_matrices()
    # inicializo la estructura
    fcc = CB.Cubic(1.465, nx=4, ny=4, nz=4, tipo="fcc")
    fcc.set_A2B2()
    # parámetros de la corrida
    corrida = PC.Corrida(1.0, 0.0001, 10000, 100, 0.5)

    fcc.write_vtf("inicio.vtf") 
    lj = LJ.LJ(components, fcc)

    energy_i = lj.energy()

    start_time = time.time()

    KT = 0.8 
    Nmc = 10000
    counter_aceptados = 0
    counter = 0
    counter_rechazados_no_iguales = 0
    counter_iguales = 0
    print counter, energy_i
    while (counter_aceptados < Nmc):
        indxs = np.random.randint(0, high=fcc.Natoms, size=2)
        a1 = fcc.IndiceElementos[indxs[0]]
        a2 = fcc.IndiceElementos[indxs[1]]
        if a1 != a2:
	    counter_iguales += 1
            if lj.distance_true(indxs[0], indxs[1]):
                deltaE = lj.energy_IJ(indxs[0], indxs[1])
                if(deltaE < 0):
                    #print (indxs, a1, a2, deltaE )
                    counter_aceptados = counter_aceptados + 1
		    energy_i +=deltaE
		    print counter, energy_i, float(counter_aceptados)/float(counter_iguales)
                else:
                    if np.exp(-KT*deltaE) > np.random.random():
                       counter_aceptados = counter_aceptados + 1
                       #print ("aceptado por KT", counter, aceptados)
		       energy_i +=deltaE
		       print counter, energy_i, float(counter_aceptados)/float(counter_iguales)
                    else:
                       pass
                       #no aceptado vuelven los índices a sus posiciones
                       fcc.IndiceElementos[indxs[0]] = a1
                       fcc.IndiceElementos[indxs[1]] = a2
        counter = counter + 1
    fcc.write_vtf("fin.vtf")
    end_time = time.time()
    print energy_i, counter, counter_aceptados, end_time - start_time


if __name__ == "__main__":
    main()
