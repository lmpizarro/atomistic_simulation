#!/usr/bin/python
# -*- coding: utf-8 -*-

#from matplotlib import mpl, pyplot
from matplotlib import pyplot
import matplotlib as mpl
import numpy as np
import os

from pylab import *
ion()

direc = 'anima'
# Abro los datos para construir la animación
with open('anim_datos.dat') as arch:
    lineas = [line.rstrip('\n') for line in arch]
    
N_tot = int(lineas[0])
mag = [float(num) for num in lineas[1:]]

archivos = [direc+'/estado%05i'%j for j in range(1,N_tot+1)]

print(np.size(archivos),np.size(mag))

# Para animacion

files=[]
filename='pelicula'

step = 0
for nom in archivos:
    data = []
    with open(nom) as ins:
        data = [[int(n) for n in line.split()] for line in ins]
    data = np.asarray(data)
    data = data[1:-1,1:-1]
    cmap = mpl.colors.ListedColormap(['red','black'])
    if step == 0:
        p = pyplot.imshow(data,interpolation='none',cmap=cmap,aspect='equal')
    else:
        p.set_data(data)
    pyplot.title('Paso:%05i'%step + '  ' + u'Magnetización:%6.1f'%mag[step])
    fname = direc + '/_anim_%05i.png'%step    
    savefig(fname)
    files.append(fname)
    step += 1
    #pyplot.pause(0.001)   
    #pyplot.show()
# Para la animaci´on
os.system("mencoder 'mf://anima/_anim_*.png' -mf type=png:fps=25 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o " + filename + ".mpg")
# cleanup
#for fname in files: os.remove(fname)