# -*- coding: utf-8 -*-

#!/usr/bin/python

import os
import subprocess
import shutil

# Número total de estados distintos que quiero guardar
N_tot = 10000
# Máximo numero de estados iguales aceptables (el programa se corta)
N_tol = 500
# Temperatura de trabajo
Temperatura = '5.5'
# Directorio donde se guardarán los archivos de los estados
direc = 'anima'

# Se crea el directorio (da error si ya existe)
os.makedirs(direc)
# Se ejecutará el programa ising pasándole datos por línea de comandos
programa = [os.getcwd() + '/ising']
# Parámetros de entrada para ising. Un paso a la temperatura definida
argumentos = ['-k'] + ['1'] + ['-t'] + [Temperatura]

j=0                 # Contador de los estados aceptados (distintos)
k=0                 # Contador de todos los estados
aux = 0             # Cantidad de estados consecutivos iguales
mag = []            # Lista con la magnetización de cada estado aceptado
eng = []            # Lista con lo energía de cada estado aceptado
mag_anterior = []   # Control para saber si cambió el estado
ind_acep = []       # Lista donde se guardan los índices de los estados aceptados

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
    # Guarda los datos de energía y magnetización
    eng.append(cols[1])
    mag.append(cols[3])
    # Se fija si cambió la magnetización
    if (cols[3] != mag_anterior) and (mag_anterior!=[]):
        # Si se produjo un cambio en la magnetización hay un nuevo estado
        print('Magnetización: %6.1f'%float(cols[3]))
        # Guarda el índice para identificar al estado diferente        
        ind_acep.append(k)
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
    k += 1
else:
    print('Se alcanzaron los {0} estados diferentes'.format(N_tot))
    
# Guardo todos los datos, por si los necesito para algo. No creo. Pero bueno.
with open('anim_datos_todos.dat','w+') as arch:
    # Total de estados    
    arch.write(str(k)+'\n')
    # Magnetización
    for num in mag:
        arch.write(str(num) +'\n')
    # Energía
    for num in eng:
        arch.write(str(num) +'\n')

# Guardo los datps que serán leidos por anima02_armado.py para hacer los
# gráficos y animaciones     
with open('anim_datos.dat','w+') as arch:
    # Total de estados    
    arch.write(str(N_tot)+'\n')
    # Magnetización
    for j in ind_acep:
        arch.write(str(mag[j]) +'\n')
    # Energía
    for j in ind_acep:
        arch.write(str(eng[j]) +'\n')
