#!/usr/bin/python
# -*- coding: utf-8 -*-


#
# calcula la fft de una archivo de velocidades
#    vx vy vz
#
import numpy as np
import matplotlib.pyplot as plt
import sys
import pyfftw as fft

arch = sys.argv[1]

line = np.loadtxt(arch,skiprows=1)
plt.figure(1,figsize=(8,6))


# determina si el archivo tiene 3 columnas
if (line.shape[1] != 3):
	sys.exit(1)

vx = np.asarray(line[:,0])
vy = np.asarray(line[:,1])
vz = np.asarray(line[:,2])

N = np.size(vx)

fvx = np.zeros(N)
fvy = np.zeros(N)
fvz = np.zeros(N)

print  "Tamaño de la muestra: ", N

# inicializa los objetos numericos para la fft 
a = fft.n_byte_align_empty(N, 16, 'complex128')
b = fft.n_byte_align_empty(N, 16, 'complex128')
c = np.zeros(N)

fft_object = fft.FFTW(a, b)

a[:] = vx + 1j*0.0
fft_a = fft_object()
fvx= np.absolute(fft_a)

a[:] = vy + 1j*0.0
fft_a = fft_object()
fvy= np.absolute(fft_a)

a[:] = vz + 1j*0.0
fft_a = fft_object()
fvz= np.absolute(fft_a)


plt.plot(fvz[1:1000])

plt.title('Archivo leido: ' + arch)
#plt.xlabel('$Tiempo$')
#plt.ylabel(u'$Presión$')
#plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)

plt.show()

print fvy, fvx, fvz
