#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script para graficar los datos del archivo tabla_temperatura.dat (o alguno con
# formato similar. Se deben cambiar los labels y títulos de acuerdo a lo que
# se quiera hacer

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    archivo = ['tabla_temperatura.dat']
else:
    archivo = sys.argv[1:]
    
forma = []
forma.append('b.')

leg = []
leg.append(ur'$\rho^* = 0.3$')

xLim = [0.65, 1.45]

j = 0
for nombre in archivo:

    data = np.loadtxt(nombre,skiprows=1)
    
    T     = data[:,0]
    rho   = data[:,1]
    p_avg = data[:,2]
    p_std = data[:,3]
    c_avg = data[:,4]
    c_std = data[:,5]
    T_avg = data[:,6]
    T_std = data[:,7]
    U_avg = data[:,8]
    U_std = data[:,9]
    
    plt.figure(1,figsize=(8,6))
    
    plt.errorbar(T_avg, p_avg, xerr=T_std, yerr=p_std, fmt=forma[j],label=leg[j])
    #plt.title(r'Valores medios')
    plt.xlabel('$T^*$')
    plt.ylabel('$p^*$')
    plt.xlim(xLim)
    plt.ylim([-0.15 , 0.25])
    plt.legend(loc='upper left',numpoints=1)
    plt.grid(True)
   
    plt.figure(2,figsize=(8,6))
    
    plt.errorbar(T_avg, c_avg, xerr=T_std,yerr=c_std, fmt= forma[j],label=leg[j])
    plt.xlabel('$T^*$')
    plt.ylabel('$\sigma(p^*)$')
    plt.xlim(xLim)
    #plt.ylim([-0.1 , 1.1])
    plt.legend(loc='upper center',numpoints=1)    
    plt.grid(True)
    
    plt.figure(3,figsize=(8,6))
    
    plt.errorbar(T_avg, U_avg, xerr=T_std,yerr=U_std, fmt= forma[j],label=leg[j])
    plt.xlabel('$T^*$')
    plt.ylabel('$U^*$')
    plt.xlim(xLim)
    plt.ylim([-4.1 , -1.9])
    plt.legend(loc='upper left',numpoints=1)    
    plt.grid(True)
#    plt.figure(2,figsize=(8,8))
#    
#    plt.subplot(2, 1, 1)
#    plt.errorbar(T, c_avg/Nsp[j], yerr=c_std/Nsp[j], fmt=forma[j],label=leg[j])
#    #plt.title('Fluctuaciones')
#    plt.xlabel('$kT/J$')
#    plt.ylabel('$<C>/N$')
#    plt.ylim([-0.1 , 2.3])
#    plt.legend(loc='upper right',numpoints=1)  
#    plt.grid(True)
#    
#    plt.subplot(2, 1, 2)
#    plt.errorbar(T, s_avg/Nsp[j], yerr=s_std/Nsp[j], fmt=forma[j],label=leg[j])
#    plt.xlabel('$kT/J$')
#    plt.ylabel('$<\chi>/N$')
#    plt.ylim([-4 , 100])
#    plt.legend(loc='upper right',numpoints=1)  
#    plt.grid(True)
    
    j += 1

plt.show()
