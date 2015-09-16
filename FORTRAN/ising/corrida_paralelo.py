#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Lo corro con:  mpirun.mpich -n 4 python corrida_paralelo.py

import os
import shutil
import subprocess
import errno
import random
import numpy as np
# Paralelización
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

# Defino la función para escribir el archivo de datos
def escribe_entrada(nombre,valor):
    with open('parametros.dat','r') as f1:
        fila = f1.readline()        
        col  = fila.split()
        if nombre== 'T':
            new=fila.replace(col[2],valor)
        elif nombre == 'N':
            new=fila.replace(col[4],valor)
        else:
            print('No se reconoce el valor a escribir')
            exit
    with open('parametros.dat','w') as fileout:
        fileout.write(new)

def escribe_semilla():
    random.seed()
    # Obtengo un número aleatorio entero de 4 bytes
    entero_aleatorio = random.randint(-2147483648,2147483647)
    with open('seed.dat','w') as fileout:
        fileout.write(str(entero_aleatorio))
        
def copia_val_medios():
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('val_medios.dat','r') as vmed:
        valor=vmed.readline()
    arch_comun = os.path.join(curr_dir,'tablas_temperatura.dat')
    with open(arch_comun,'a') as comun:
        comun.write(valor)        

# Información de dónde está            
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()
# Lista con las temperaturas que se desean calcular
T_min = 0.5
T_max = 1.2
dT    = 0.1
tempe = np.arange(T_min,T_max+dT,dT)
#tempe = tempe.tolist() + [2.15, 2.25, 2.35, 2.45, 2.55] # Mäs detalle en la Tc
tempe.sort()
# lo convierto a string
tempe = [str(i) for i in tempe]

# Pralelización - Separo lista en la cantidad de cores que quiero correr
nT_local = len(tempe)/size      # Debe ser entero

tempe_local=tempe[rank*nT_local:(rank+1)*nT_local]

print('El core',rank, 'hace las corrida de las temperaturas', tempe_local)

# Loop para crear todos los directorios y correr el ejecutable en ellos
for T in tempe_local:
    # Nombre de la carpeta uqe se va a crear
    carpeta = T + '_tmpfolder'
    # Camino completo de la carpeta que se va a crear
    path_carpeta = os.path.join(curr_dir,carpeta)
    # Se crea la carpeta sólo si ésta no existe. De lo contrario se saltea la 
    # temperatura correspondiente
    try: 
        os.makedirs(path_carpeta)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print('Ya existe el directorio. Se omite corrida con T=' + T)
            continue
    # Copia el archivo de entrada a la carpeta 
    shutil.copy('parametros.dat',path_carpeta)
    #shutil.copy('ising',path_carpeta) 
    # Se mete en la carpeta
    os.chdir(path_carpeta)       
    # Cambia el archivo de entrada adentro de la carpeta
    escribe_entrada('T',T)
    # Escribe la semmilla en la carpeta
    escribe_semilla()

    # Corre el programa para ver la convergencia
    escribe_entrada('N','200000') 
    salida = subprocess.check_output(curr_dir+'/ising')
    
    # Guardo la salida para ver ue hizo
    f=open("log1.txt",'w')
    f.write(salida)
    f.close() 


#############################
#initializing variables. mpi4py requires that we pass numpy objects.
aviso = np.ones(1)
listos_buffer = np.zeros(1)


# root node receives results from all processes and sums them
if rank == 0:
    listos = listos_buffer[0]
    for i in range(1, size):
        comm.Recv(listos_buffer, ANY_SOURCE)
        listos += listos_buffer[0]
else:
        # all other process send their result
        comm.Send(aviso)
        
if comm.rank == 0:
        print(listos)
        
 # Mando la segunda corrida a cada directorio       
for T in tempe_local:
    # Nombre de la carpeta uqe se va a crear
    carpeta = T + '_tmpfolder'
    # Camino completo de la carpeta que se va a crear
    path_carpeta = os.path.join(curr_dir,carpeta)

    # Se mete en la carpeta
    os.chdir(path_carpeta)       
    # Cambia el archivo de entrada adentro de la carpeta
    escribe_entrada('T',T)
    # Escribe la semmilla en la carpeta
    escribe_semilla()

    # Corre por segunda vez tomando el estado anterior. Aumento el N
    escribe_entrada('N','5000000')
    salida = subprocess.check_output(curr_dir+'/ising')
    
    # Guardo la salida para ver ue hizo
    f=open("log2.txt",'w')
    f.write(salida)
    f.close()    

    # Copia los valores medios
    copia_val_medios()
    # Sale de la carpeta
    os.chdir(curr_dir)