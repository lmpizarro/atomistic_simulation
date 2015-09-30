# -*- coding: utf-8 -*-

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch_eng = 'energia.dat'
    arch_mag = 'magneti.dat'
else:
    arch_eng = sys.argv[1]
    arch_mag = sys.argv[1]


E = np.loadtxt(arch_eng)
M = np.loadtxt(arch_mag)


plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(E/400, 'b.')
plt.title('Datos temporales')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Energía')

plt.subplot(2, 1, 2)
plt.plot(M/400, 'b.')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Magnetización')


plt.show()
