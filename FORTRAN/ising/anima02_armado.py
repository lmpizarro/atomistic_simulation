#!/usr/bin/python
# -*- coding: utf-8 -*-

from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter # Cambia formato labels grafico
import numpy as np
import os
import sys

from pylab import *
#ion()

# Directorio donde se guardaron los estados que se desean graficar (realizado
# con el script anima01_corridas.py)
direc = 'anima'
# Nombre del archivo de video que se grabará
filename='pelicula'

###############################################################################
# Lee argumentos por línea de comandos para saber si graba video
###############################################################################

if ('grabar' in sys.argv):
    graba = True
else:
    graba = False


###############################################################################
# Leo estados producidos por el script anima01_corridas.py
###############################################################################

# Abra los datos para construir la animación
with open('anim_datos.dat') as arch:
    lineas = [line.rstrip('\n') for line in arch]
    
# Número total de imágenes
N_tot = int(lineas[0])
# Valores de magnetización
mag = [float(num) for num in lineas[1:N_tot+1]]
# Valores de energía
eng = [float(num) for num in lineas[N_tot+1:]]
# Archivos con los estados
archivos = [direc+'/estado%05i'%j for j in range(1,N_tot+1)]

###############################################################################
# Comienza la graficación y animación
###############################################################################
files=[]    # Nombres de los archivos png (para borrarlos al final)
step = 0    # Contador
gs = gridspec.GridSpec(4, 2,
                       width_ratios=[1,1],
                       height_ratios=[20,1,5,5]
                       )
majorFormatter = FormatStrFormatter('%5.f')
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
        f1 = pyplot.figure(0,figsize=(7,8))
        # Matriz de spines
        ax0=f1.add_subplot(gs[0:-3,:])
        p1 = pyplot.imshow(data,interpolation='none',cmap=cmap,aspect='equal')
        # Energía        
        ax1 = f1.add_subplot(gs[-2,:])        
        ax1.grid(True)
        ax1.set_title(u'Energía',color='b')
        p2a, = ax1.plot(step,eng[step],'b.-')
        ax1.xaxis.set_ticklabels([])
        # Magnetización        
        ax2 = f1.add_subplot(gs[-1,:]) 
        ax2.set_xlabel('Paso')
        ax2.grid(True)
        ax2.set_title(u'Magnetización',color='g')
        p2b, = ax2.plot(step,mag[step],'g.-')
    else:
        # Actualizo la figura
        p1.set_data(data)
        # Actualiza energía
        p2a.set_data(range(step),eng[0:step])
        ax1.relim()
        ax1.autoscale_view(True,True,True)
        ax1.yaxis.set_major_formatter(majorFormatter)
        # Actualiza magnetización        
        p2b.set_data(range(step),mag[0:step])
        ax2.relim()
        ax2.autoscale_view(True,True,True)    
        ax2.yaxis.set_major_formatter(majorFormatter)
        # Dibuja (en el caso de graficar todo junto)        
        pyplot.draw()

    # Título con información relevante
    titulo = 'Paso:%05i'%step + '   ' + \
             'Mag:%6.1f'%mag[step] + '   '+ \
             'Eng:%6.1f'%eng[step]
    ax0.set_title(titulo)
    # Nombre de los archivos png que voy a guardar para hacer un video    
    if graba:    
        fname = direc + '/_anim_%05i.png'%step    
        savefig(fname)
        files.append(fname)
    step += 1
    pyplot.pause(0.00001)   
    #pyplot.show()

###############################################################################
# Para construir el video se utiliza el programa externo mencoder
###############################################################################
# Opciones que se le pasan a mencoder
if graba:
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