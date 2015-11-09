#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch = 'gr.dat'
else:
    arch_tem = sys.argv[1]

# Se leen los dos archivos, salvo la primer linea
datos = np.loadtxt(arch,skiprows=1)

r  = datos[:,0]
gr = datos[:,1]

plt.figure(1,figsize=(8,6))

plt.plot(r, gr, 'k.')
#plt.title('Datos temporales')
plt.xlabel(u'$r^*$')
plt.ylabel(u'$g(r^*)$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)


plt.grid(True)

plt.show()
