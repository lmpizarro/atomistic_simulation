#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch = 'gr.dat'
else:
    arch = sys.argv[1:]

forma = ['g.', 'b.', 'r.', 'k.']
#forma.append('b.')

# if len(arch)==2:
#     leg = [r'$T^*=1.1 \quad \rho^* = 0.3$',
#            r'$T^*=1.4 \quad \rho^* = 0.3$']
# elif len(arch)==3:
#     leg = [r'$T^*=1.1 \quad \rho^* = 0.001$',
#            r'$T^*=1.1 \quad \rho^* = 0.3$',
#            r'$T^*=1.1 \quad \rho^* = 0.8$']

leg = [r'$g_{11}(r)$',
       r'$g_{12}(r)$',
       r'$g_{22}(r)$',
       r'$g(r)$']

#xlim = [0.0 , 4.0]

j = 0

#for nombre in arch:
    # Se leen los dos archivos, salvo la primer linea
datos = np.loadtxt(arch,skiprows=1)

r   = datos[:,0]
d   = datos.shape
gr  = np.zeros((d[0],d[1]-1))

for j in range(4):
    gr[:,j] = datos[:,j+1]
    plt.figure(1,figsize=(8,6))
    plt.plot(r, gr[:,j], forma[j], label=leg[j])
     #   plt.xlim(xlim)
    plt.grid(True)
    plt.legend(loc='upper right',numpoints=1)

plt.xlabel(u'$r^*$')
plt.ylabel(u'$g(r^*)$')

plt.show()
