# -*- coding: utf-8 -*-

#!/usr/bin/python

# Scripts para graficar los archivos con las correlaciones producidas
# con el otro script <correlaciones.py>. 
# Se le debe psar el nombre de los archivos con los datos
# Está armado paara que se le ingresen tres archivos de distinta
# temperatura de una red de 20x20. De lo contrario hay que cambiar
# labels y títulos.

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    archivo = ['energia.dat']
else:
    archivo = sys.argv[1:]


forma = []
forma.append('k.')
forma.append('b.')
forma.append('g.')

leg = []
leg.append(r'$T = 1.0$')
leg.append(r'$T = 2.3$')
leg.append(r'$T = 9.0$')


j = 0
for nombre in archivo:

    data = np.loadtxt(nombre)

    t = data[:,0]
    c_avg = data[:,1]
    c_std = data[:,2]

    plt.figure(1,figsize=(8,8))

    # Zona de  graficaciones log
    n_log = 400

    plt.subplot(2, 1, 1)
    #plt.plot(t,corr_avg, 'b.')
    plt.errorbar(t, c_avg, yerr=c_std,fmt=forma[j],label=leg[j])
    plt.xlabel('$k-Paso\, MC$')
    plt.ylabel(u'$C_{EE}[k]$')
    plt.legend(loc='upper center',numpoints=1)
    plt.title('Autocorrelaciones para red de 20 x 20')
    plt.grid(True) 
    # Zona de ajuste
    n_fin = 250

    p = np.polyfit(t[:n_fin],np.log(c_avg[:n_fin]),1)
    t_corr =- 1/p[0]
    plt.subplot(2,1,2)
    corr_E_aj = np.poly1d(p)
    plt.plot(t[:n_fin],corr_E_aj(t[:n_fin]),'r-')
    plt.ylim(-2.0,0.1)
    plt.grid(True)


    plt.subplot(2, 1, 2)
    plt.plot(t[:n_log],np.log(c_avg[:n_log]), forma[j],label=r'$\tau_c = {:.1f}$'.format(t_corr))
    plt.xlabel('$k-Paso\, MC$')
    plt.ylabel(u'$ln(C_{EE}[k])$')
    #plt.legend(loc='upper right',numpoints=1)
    plt.legend(loc='lower left',numpoints=1)

    plt.plot(t[:n_fin],corr_E_aj(t[:n_fin]),'r-')
    j += 1
   
plt.show()
