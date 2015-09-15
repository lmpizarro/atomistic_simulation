# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

archivo = 'tablas_temperatura.dat'


data = np.loadtxt(archivo)

T = data[:,0]
E_avg = data[:,1]
E_var = data[:,2]
M_avg = abs(data[:,3])
M_var = data[:,4]

plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(T, E_avg/400, 'k.')
plt.title('Temperatura')
plt.ylabel('Energia media')

plt.subplot(2, 1, 2)
plt.plot(T, M_avg/400, 'k.')
plt.xlabel('Temperatura')
plt.ylabel('Magnetizacion media')

plt.figure(2)

plt.subplot(2, 1, 1)
plt.plot(T, E_var, 'k.')
plt.title('Temperatura')
plt.ylabel('Varianza de la energia')

plt.subplot(2, 1, 2)
plt.plot(T, M_var, 'k.')
plt.xlabel('Temperatura')
plt.ylabel('Varianza de la magnetizacion')

plt.show()
