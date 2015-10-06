# -*- coding: utf-8 -*-

#!/usr/bin/python

# Script para calcular la autocorrelación
# se le debe ingresar el archivo con los valores de energía
# Guarda la salida en otro archivo para ser leido luego
# por los scripts de graficación

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

n_por_block = 6000
N_block     = N_E // n_por_block


# Si no quiero hacer muchos promedios, tomo una cantidad de bloques fija
N_minima = 100
N_block = min(N_minima,N_block)

E_acum=0

# Define la matriz donde se guardarán las correlaciones
corr_tot = np.zeros((N_block,n_por_block))


for j in range(N_block):
    E_block = E[j*n_por_block:(j+1)*n_por_block]
    E_block = E_block - np.mean(E_block)
    # Molesto si tengo valores constantes
    if np.std(E_block)==0:
        norm_E = 1
        E_block = np.ones(np.size(E_block))
    else:
        norm_E = np.std(E_block)**2
    print('Calculando correlación del bloque {0}'.format(j))
    corr_E_block = signal.correlate(E_block,E_block,mode='full')
    corr_tot[j][:] = corr_E_block[n_por_block-1:]/norm_E

nor_bias = np.array([n_por_block-j for j in range(n_por_block)])

corr_tot = corr_tot/nor_bias
corr_avg = np.mean(corr_tot,0)
corr_std = np.std(corr_tot,0)/np.sqrt(N_block)

t = np.arange(np.size(corr_avg))

plt.figure(1)

# Zona de  graficaciones log
n_log = 400

plt.subplot(2, 1, 1)
plt.plot(t[:n_log],np.log(corr_avg[:n_log]), 'b.')
plt.title('Datos temporales')
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Correlacion Energía')

plt.subplot(2, 1, 2)
#plt.plot(t,corr_avg, 'b.')
plt.errorbar(t, corr_avg, yerr=corr_std)
plt.xlabel('Tiempo [u.a.]')
plt.ylabel(u'Correlacion Energia')
#plt.yscale('log')

# Zona de ajuste
n_fin = 200

p = np.polyfit(t[:n_fin],np.log(corr_avg[:n_fin]),1)
t_corr =-1/p[0]
plt.subplot(2,1,1)
corr_E_aj = np.poly1d(p)
plt.plot(t[:n_fin],corr_E_aj(t[:n_fin]),'r-')
#plt.xlim(0,2*t[n_fin])
plt.title(u'Tiempo de correlación = {}'.format(t_corr))

# Para escribir en un archivo
lis = [t, corr_avg, corr_std]

archi_escr = 'corre_'+ arch_eng

with open(archi_escr,'w') as f:
    for x in zip(*lis):
        f.write('{0}\t{1}\t{2}\n'.format(*x))

plt.show()
