# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

arch_eng = 'energia.dat'
arch_mag = 'magneti.dat'

E = np.loadtxt(arch_eng)
M = np.loadtxt(arch_mag)


plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(E/400, 'b.')
plt.title('Datos temporales')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel('Energia')

plt.subplot(2, 1, 2)
plt.plot(M/400, 'b.')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel('Magnetizacion')


plt.show()
