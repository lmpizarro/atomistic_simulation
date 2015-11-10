#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

plt.figure(1,figsize=(8,6))

arch = sys.argv[1]

y = np.loadtxt(arch,skiprows=1)

if len(y.shape) > 1:
    columnas = y.shape[1]
else:
    columnas = 1

if columnas == 2:
    x = y[:,0]
    y = y[:,1] 
    plt.plot(x,y, 'k.')
else:
    plt.plot(y, 'k.')


plt.title('Archivo leido: ' + arch)
#plt.xlabel('$Tiempo$')
#plt.ylabel(u'$Presión$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)

plt.show()
