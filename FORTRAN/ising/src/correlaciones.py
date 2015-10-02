# -*- coding: utf-8 -*-

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import signal

if len(sys.argv) < 2:
    arch_eng = 'energia.dat'
    arch_mag = 'magneti.dat'
else:
    arch_eng = sys.argv[1]
    arch_mag = sys.argv[1]


E = np.loadtxt(arch_mag)
#M = np.loadtxt(arch_mag)

N_E = np.size(E)
#N_M = np.size(M)

n_por_block = 5000
N_block     = N_E // n_por_block


# Si no quiero hacer muchos promedios, tomo una cantidad de bloques fija
N_minima = 100
N_block = min(N_minima,N_block)

E_acum=0

for j in range(N_block):
    E_block = E[j*n_por_block:(j+1)*n_por_block]
    E_block = E_block - np.mean(E_block)
    norm_E = np.std(E_block)**2
    print('Calculando correlación del bloque {0}'.format(j))
    corr_E_block = signal.correlate(E_block,E_block,mode='full')
    E_acum = E_acum + corr_E_block/norm_E

print(np.size(E_acum))

corr_E = E_acum[n_por_block:]/N_minima/n_por_block


t = np.arange(np.size(corr_E))

plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(t,np.log(corr_E), 'b.')
plt.title('Datos temporales')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Correlacion Energía')

plt.subplot(2, 1, 2)
plt.plot(t,corr_E, 'b.')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Correlacion Energia')
#plt.yscale('log')

# Zona de ajuste
n_fin = 200

p = np.polyfit(t[:n_fin],np.log(corr_E[:n_fin]),1)
t_corr = 1/p[1]
plt.subplot(2,1,1)
corr_E_aj = np.poly1d(p)
plt.plot(t[:n_fin],corr_E_aj(t[:n_fin]),'r-')
#plt.xlim(0,2*t[n_fin])
plt.title(u'Tiempo de correlación = {}'.format(t_corr))


plt.show()
