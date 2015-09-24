#!/usr/bin/python
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
E_std = data[:,2]
M_avg = np.abs(data[:,3])
M_std = data[:,4]
c_avg = data[:,5]
c_std = data[:,6]
s_avg = data[:,7]
s_std = data[:,8]
a_avg = data[:,9]
a_std = data[:,10]

plt.figure(1)

plt.subplot(2, 1, 1)
plt.errorbar(T, E_avg/400, yerr=E_std/400, fmt='k.')
plt.title('Valores medios')
plt.xlabel('Temperatura')
plt.ylabel(u'Energía')
plt.ylim([-2.1 , -0.2])

plt.subplot(2, 1, 2)
plt.errorbar(T, M_avg/400,yerr=M_std/400, fmt= 'k.')
#M_teo = ( 1 - (np.sinh(2/T))**(-4) )**(1/8)
M_teo = ( 1 - (np.sinh(2/T))**(-4) )**(float(1)/8)
plt.plot(T,M_teo,'r-')
plt.xlabel('Temperatura')
plt.ylabel(u'Magnetización')
plt.ylim([-0.1 , 1.1])
#plt.ylim([-1.1 , 0.1])

plt.figure(2)

plt.subplot(2, 1, 1)
plt.errorbar(T, c_avg/400, yerr=c_std/400, fmt='k.')
plt.title('Fluctuaciones')
plt.xlabel('Temperatura')
plt.ylabel(u'Calor específico')
plt.ylim([-0.1 , 1.9])

plt.subplot(2, 1, 2)
plt.errorbar(T, s_avg/400, yerr=s_std/400, fmt='k.')
plt.xlabel('Temperatura')
plt.ylabel(u'Suceptibilidad magnética')
plt.ylim([-4 , 60])

plt.figure(3)

plt.errorbar(T,a_avg,yerr=a_std,fmt='k.')
plt.xlabel('Temperatura')
plt.ylabel('Aceptaciones [%]')
plt.ylim([-5 , 75])

plt.show()
