# -*- coding: utf-8 -*-

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    arch_eng = 'energia.dat'
    arch_mag = 'magneti.dat'
else:
    arch_eng = sys.argv[1]
    arch_mag = sys.argv[1]


E = np.loadtxt(arch_eng)
M = np.loadtxt(arch_mag)


plt.figure(1,figsize=(8,8))

plt.subplot(2, 1, 1)
plt.plot(E/400, 'k.')
#plt.title('Datos temporales')
plt.xlabel('$Paso\,MC$')
plt.ylabel(u'$E/N$')
plt.ylim(-2.1,plt.ylim()[1])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(M/400, 'k.')
plt.xlabel('$Paso\,MC$')
plt.ylabel(u'$M/N$')
if M[-1]/400 < 0:
    plt.ylim(-1.1,plt.ylim()[1])
else:
    plt.ylim(plt.ylim()[0],1.1)

plt.grid(True)

plt.show()
