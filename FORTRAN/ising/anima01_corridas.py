#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess
import shutil

# Número total de estados distintos que quiero guardar
N_tot = 10000
# Máximo numero de estados iguales aceptables (el programa se corta)
N_tol = 500
# Temperatura de trabajo
Temperatura = '0.7'
# Directorio donde se guardarán los archivos de los estados
direc = 'anima'

# Se crea el directorio (da error si ya existe)
os.makedirs(direc)
# Se ejecutará el programa ising pasándole datos por línea de comandos
programa = [os.getcwd() + '/ising']
# Parámetros de entrada para ising. Un paso a la temperatura definida
argumentos = ['-k'] + ['1'] + ['-t'] + [Temperatura]

j=0                 # Contador de los estados aceptados (distintos)
aux = 0             # Cantidad de estados consecutivos iguales
mag = []            # Lista con la magnetización de cada estado aceptado
mag_anterior = []   # Control para saber si cambió el estado

# Se ejecuta hasta alcanzar los N_tot, o hasta que pasan aux estados sin cambio

while (j<N_tot):
    # Ejecuta ising
    proc = subprocess.Popen(programa + argumentos,stdout=subprocess.PIPE)
    # Espera que termine    
    proc.wait()
    # Lee el estado en que quedó del archivo val_medios.dat
    with open('val_medios.dat') as arch:
        line = arch.readline()
        cols = line.split()
    # Se fija si cambió la magnetización
    if (cols[3] != mag_anterior) and (mag_anterior!=[]):
        # Si se produjo un cambio en la magnetización hay un nuevo estado
        print('Magnetización: %6.1f'%float(cols[3]))
        # Guardo la magnetización        
        mag.append(cols[3])        
        j += 1
        # Copio el estado final
        shutil.copy('estado_mat.dat','./'+direc+'/estado%05i'%j)
        # Hago que el estado final sea el nuevo estado inicial        
        os.rename('estado.dat','ultimo_estado.dat')
        # Guardo magnetización para comparar
        mag_anterior = cols[3]
        # Reinicializo el contado
        aux = 0
    else:
        # Si el estado no cambió
        mag_anterior = cols[3]
        # Cuento la cantidad de estados sin cambiar
        aux += 1
        # Si llego a un número igual a N_tol, se corta el programa
        if aux > N_tol:
            print('No cambia luego de {0} iteraciones'.format(N_tol))
            N_tot = j
            break
else:
    print('Se alcanzaron los {0} estados diferentes'.format(N_tot))
    
# Guardo los datps que serán leidos por anima02_armado.py para hacer los
# gráficos y animaciones
with open('anim_datos.dat','w+') as arch:
    arch.write(str(N_tot)+'\n')
    for num in mag:
        arch.write(str(num) +'\n')
