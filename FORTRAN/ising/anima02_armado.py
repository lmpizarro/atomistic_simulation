#!/usr/bin/python
# -*- coding: utf-8 -*-

from matplotlib import pyplot
import matplotlib as mpl
import numpy as np
import os

from pylab import *
ion()


# Directorio donde se guardaron los estados que se desean graficar (realizado
# con el script anima01_corridas.py)
direc = 'anima'
# Nombre del archivo de video que se grabará
filename='pelicula'

##############################################################################
# Leo estados producidos por el script anima01_corridas.py
##############################################################################

# Abra los datos para construir la animación
with open('anim_datos.dat') as arch:
    lineas = [line.rstrip('\n') for line in arch]
    
# Número total de imágenes
N_tot = int(lineas[0])
# Valores de magnetización
mag = [float(num) for num in lineas[1:]]
# Archivos con los estados
archivos = [direc+'/estado%05i'%j for j in range(1,N_tot+1)]

##############################################################################
# Comienza la graficación y animación
##############################################################################
files=[]    # Nombres de los archivos png (para borrarlos al final)
step = 0    # Contador

for nom in archivos:
    data = []
    # Abre el archivo con los estados que se quieren graficar
    with open(nom) as ins:
        data = [[int(n) for n in line.split()] for line in ins]
    data = np.asarray(data)    # Lo convierte a array de numpy para usar slices
    data = data[1:-1,1:-1]     # Saco los bordes (que se puesieron pro las CC)
    # Mapeo de los -1 y 1 a colores    
    cmap = mpl.colors.ListedColormap(['red','black']) 
    if step == 0:
        # Sólo en el primer paso para hacerlo más rápido
        p = pyplot.imshow(data,interpolation='none',cmap=cmap,aspect='equal')
    else:
        # Actualizo la figura
        p.set_data(data)
    # Título con información relevante
    pyplot.title('Paso:%05i'%step + '  ' + u'Magnetización:%6.1f'%mag[step])
    # Nombre de los archivos png que voy a guardar para hacer un video    
    fname = direc + '/_anim_%05i.png'%step    
    savefig(fname)
    files.append(fname)
    step += 1
    pyplot.pause(0.0001)   
    #pyplot.show()

##############################################################################
# Para construir el video se utiliza el programa externo mencoder
##############################################################################
# Opciones que se le pasan a mencoder
cmdstring = ('mencoder ' +
            'mf://%s/_anim_*.png '%direc +
            '-mf ' + 'type=png:fps=40 '+
            '-ovc '+ 'lavc ' + '-lavcopts '+
            'vcodec=wmv2 '+
            '-oac ' + 'copy '+
            '-o '+filename+'.mpg')

os.system(cmdstring)

# Una opción para comprimir el video es:
#mencoder -oac mp3lame -ovc x264 pelicula.mpg -o pelicula.avi
# Borra todos los archivos png
#for fname in files: os.remove(fname)