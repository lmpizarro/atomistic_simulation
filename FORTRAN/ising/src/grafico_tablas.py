# -*- coding: utf-8 -*-

#!/usr/bin/python

# Script para graficar los datos del archivo tablas_temperatura (o alguno con
# formato similar. Se deben cambiar los labels y títulos de acuerdo a lo que
# se quiera hacer

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    archivo = ['tablas_temperatura.dat']
else:
    archivo = sys.argv[1:]
    
Nsp = []
Nsp.append(20*20)

forma = []
forma.append('k.')

leg = []
leg.append('20 x 20')

j = 0
for nombre in archivo:

    data = np.loadtxt(nombre)
    
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
    
    plt.figure(1,figsize=(8,8))
    
    plt.subplot(2, 1, 1)
    plt.errorbar(T, E_avg/Nsp[j], yerr=E_std/Nsp[j], fmt=forma[j],label=leg[j])
    #plt.title(r'Valores medios')
    plt.xlabel('$kT/J$')
    plt.ylabel('$<E/J>/N$')
    plt.ylim([-2.1 , -0.2])
    plt.legend(loc='upper left',numpoints=1)
    plt.grid(True)
   
    
    plt.subplot(2, 1, 2)
    plt.errorbar(T, M_avg/Nsp[j],yerr=M_std/Nsp[j], fmt= forma[j],label=leg[j])
    plt.xlabel('$kT/J$')
    plt.ylabel('$<M/\mu>/N$')
    plt.ylim([-0.1 , 1.1])
    plt.legend(loc='upper right',numpoints=1)    
    plt.grid(True)
    
    plt.figure(2,figsize=(8,8))
    
    plt.subplot(2, 1, 1)
    plt.errorbar(T, c_avg/Nsp[j], yerr=c_std/Nsp[j], fmt=forma[j],label=leg[j])
    #plt.title('Fluctuaciones')
    plt.xlabel('$kT/J$')
    plt.ylabel('$<C>/N$')
    plt.ylim([-0.1 , 2.3])
    plt.legend(loc='upper right',numpoints=1)  
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.errorbar(T, s_avg/Nsp[j], yerr=s_std/Nsp[j], fmt=forma[j],label=leg[j])
    plt.xlabel('$kT/J$')
    plt.ylabel('$<\chi>/N$')
    plt.ylim([-4 , 100])
    plt.legend(loc='upper right',numpoints=1)  
    plt.grid(True)
    
    plt.figure(3)
    
    plt.errorbar(T,a_avg,yerr=a_std,fmt=forma[j],label=leg[j])
    plt.xlabel('$kT/J$')
    plt.ylabel('$Aceptaciones [\%]$')
    plt.ylim([-5 , 75])
    plt.legend(loc='upper left',numpoints=1)
    plt.grid(True)
    
    j += 1

plt.show()
