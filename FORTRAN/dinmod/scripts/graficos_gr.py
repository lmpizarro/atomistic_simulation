#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch = 'gr.dat'
else:
    arch = sys.argv[1:]

forma = ['k.', 'b.', 'r.', 'g.']
forma.append('b.')

if len(arch)==2:
    leg = [r'$T^*=1.1 \quad \rho^* = 0.3$',
           r'$T^*=1.4 \quad \rho^* = 0.3$']
elif len(arch)==3:
    leg = [r'$T^*=1.1 \quad \rho^* = 0.001$',
           r'$T^*=1.1 \quad \rho^* = 0.3$',
           r'$T^*=1.1 \quad \rho^* = 0.8$']

xlim = [0.0 , 4.0]

j = 0

for nombre in arch:
    # Se leen los dos archivos, salvo la primer linea
    datos = np.loadtxt(nombre,skiprows=1)

    r  = datos[:,0]
    gr = datos[:,1]

    plt.figure(1,figsize=(8,6))

    plt.plot(r, gr, forma[j], label=leg[j])
    plt.xlabel(u'$r^*$')
    plt.ylabel(u'$g(r^*)$')
    plt.xlim(xlim)
    plt.grid(True)
    plt.legend(loc='upper right',numpoints=1)


    j += 1

plt.show()
