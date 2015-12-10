#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

plt.figure(1,figsize=(8,6))

arch = sys.argv[1]

datos = np.loadtxt(arch,skiprows=1)

x = datos[:,0]
y = datos[:,1] 
z = datos[:,2]
plt.plot(x, 'k.-',label='psd(v_x)')
plt.plot(y, 'r.-',label='psd(v_y)')
plt.plot(z, 'b.-',label='psd(v_z)')

plt.title('Archivo leido: ' + arch)
plt.xlabel('$Paso$')
#plt.ylabel(u'$Presión$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)
plt.legend(loc='upper center',numpoints=1)

plt.show()
