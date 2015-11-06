#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

arch = sys.argv[1]

y = np.loadtxt(arch,skiprows=1)

plt.figure(1,figsize=(8,6))

plt.plot(y, 'k.')
plt.title('Archivo = ' + arch)
#plt.xlabel('$Tiempo$')
#plt.ylabel(u'$Presión$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)

plt.show()
