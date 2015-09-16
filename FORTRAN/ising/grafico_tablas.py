# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    archivo = 'tablas_temperatura.dat'
else:
    archivo = sys.argv[1]

data = np.loadtxt(archivo)

T = data[:,0]
E_avg = data[:,1]
E_var = data[:,2]
M_avg = abs(data[:,3])
M_var = data[:,4]

plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(T, E_avg/400, 'b.')
plt.title('Valores medios')
plt.xlabel('Temperatura')
plt.ylabel('Energia media')

plt.subplot(2, 1, 2)
plt.plot(T, M_avg/400, 'b.')
plt.xlabel('Temperatura')
plt.ylabel('Magnetizacion media')

plt.figure(2)

plt.subplot(2, 1, 1)
plt.plot(T, E_var, 'b.')
plt.title('Varianzas')
plt.xlabel('Temperatura')
plt.ylabel('Varianza de la energia')

plt.subplot(2, 1, 2)
plt.plot(T, M_var, 'b.')
plt.xlabel('Temperatura')
plt.ylabel('Varianza de la magnetizacion')

plt.show()
