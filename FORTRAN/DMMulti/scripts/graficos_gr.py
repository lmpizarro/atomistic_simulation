#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch = 'gr.dat'
else:
    arch = sys.argv[1:]

forma = ['g.-', 'b.-', 'r.-', 'k.-']
#forma.append('b.')

# if len(arch)==2:
#     leg = [r'$T^*=1.1 \quad \rho^* = 0.3$',
#            r'$T^*=1.4 \quad \rho^* = 0.3$']
# elif len(arch)==3:
#     leg = [r'$T^*=1.1 \quad \rho^* = 0.001$',
#            r'$T^*=1.1 \quad \rho^* = 0.3$',
#            r'$T^*=1.1 \quad \rho^* = 0.8$']

leg = [r'$g_{11}(r)$',
       r'$g_{12}(r)$',
       r'$g_{22}(r)$',
       r'$g(r)$']

xlim = [0.5 , 2.0]

j = 0

#for nombre in arch:
    # Se leen los dos archivos, salvo la primer linea
datos = np.loadtxt(arch,skiprows=1)

r   = datos[:,0]
d   = datos.shape
gr  = np.zeros((d[0],d[1]-1))

fig= range(4)
ax = range(4)

#for j in range(4):
#    fig[j] = plt.figure(1,figsize=(8,6))
#    ax[j]  = fig[j].add_subplot(111)    
    

fig1 = plt.figure(1,figsize=(8,6))
ax1  = fig1.add_subplot(111)
ax1.set_xlabel(u'$r^*$')
ax1.set_ylabel(u'$g(r^*)$')
ax1.grid(True)
fig2 = plt.figure(2,figsize=(8,6))
ax2  = fig2.add_subplot(111)
ax2.set_xlabel(u'$r^*$')
ax2.set_ylabel(u'$g(r^*)$')
ax2.grid(True)
fig3 = plt.figure(3,figsize=(8,6))
ax3  = fig3.add_subplot(111)
ax3.set_xlabel(u'$r^*$')
ax3.set_ylabel(u'$g(r^*)$')
ax3.grid(True)

for j in range(4):
    gr[:,j] = datos[:,j+1]
for j in [0,2]:
    ax1.plot(r, gr[:,j], forma[j], label=leg[j])
    ax1.legend(loc='upper right',numpoints=1)
    ax1.set_xlim(xlim)
for j in [1,3]:
    ax2.plot(r,gr[:,j],forma[j],label=leg[j])
    ax2.legend(loc='upper right',numpoints=1)
    ax2.set_xlim(xlim)
for j in [0,1,2,3]:
    ax3.plot(r,gr[:,j],forma[j],label=leg[j])
    ax3.legend(loc='upper right',numpoints=1)
    ax3.set_xlim(xlim)


plt.show()
