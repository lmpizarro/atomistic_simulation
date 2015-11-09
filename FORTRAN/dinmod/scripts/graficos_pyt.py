#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch_pre = 'presion.dat'
    arch_tem = 'temperatura.dat'
else:
    arch_pre = sys.argv[1]
    arch_tem = sys.argv[1]

# Se lee la primer linea
# Para la presión
with open(arch_pre) as f:
    N,dt = f.readline().split()
N_p  = np.int(N)
dt_p = np.float(dt)
# Para la temperatura
with open(arch_tem) as f:
    N,dt = f.readline().split()
N_t  = np.int(N)
dt_t = np.float(dt)

# Vector temporal
tie_p  = np.arange(0,N_p*dt_p,dt_p)
tie_t  = np.arange(0,N_t*dt_t,dt_t)
# Se leen los dos archivos, salvo la primer linea
pres = np.loadtxt(arch_pre,skiprows=1)
temp = np.loadtxt(arch_tem,skiprows=1)

plt.figure(1,figsize=(8,8))

plt.subplot(2, 1, 1)
plt.plot(tie_p, pres, 'k.')
#plt.title('Datos temporales')
plt.xlabel('$t^*$')
plt.ylabel(u'$p^*$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(tie_t, temp, 'k.')
plt.xlabel('$t^*$')
plt.ylabel(u'$T^*$')

plt.grid(True)

plt.show()
