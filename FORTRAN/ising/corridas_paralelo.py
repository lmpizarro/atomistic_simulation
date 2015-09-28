#!/usr/bin/python
# -*- coding: utf-8 -*-
 
###############################################################################       
#  SE EJECUTA CON: mpirun.mpich -n 4  python corridas_paralelo.py
###############################################################################
#
# Script para correr el programa de ising a distintas temperaturas
# Crea una carpeta para cada temperatura. En cada una de esas carpetas, a su
# vez, crea Nruns carpetas para hacer estadística y obtener los valores con
# sus respectivos errores.
#
# En cada temperatura, se utiliza el valor final de la temperatura anterior.
# Arbitrariamente, se toma el valor final de RUN00 como el valor inicial de
# todas las corridas a la siguiente temperatura (menor)
#
# Los resultados los guarda en el archivo 'tablas_temperatura.dat'
# Cada fila de 'tablas_temperatura.dat' representa los resultados a una dada
# temperatura con sus respectivos errores. Las columnas son:
# T <E> std(E) <M> std(M) <c> std(c) <suc> std(suc) <acept> std(acept)
# Los errores se calculan como std(RUN)/sqrt(NRUN)
#
# Se paralelizan las corridas que se hace a una dada tempteratura. La limitación
# está en que se debe esperar que terminen las corridas a la temperatura anterior
# para utilizar el estado final. De lo contrario se podría parelelizar también
# la corrida para cada temperatura (en el archivo 'sin errores' está hecho así)
###############################################################################

import os
import shutil
import subprocess
import errno
import numpy as np

# Funciones propias para este script
import ising_fun as isf

# Paralelización
from mpi4py import MPI
#from mpi4py.MPI import ANY_SOURCE

# Información de dónde está            
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

###############################################################################       
#   PARAMETROS DE ENTRADA
###############################################################################
# Barrido de temperaturas
# Temperatura mínima
T_min = np.float(0.5)
# Temperatura máxima
T_max = np.float(4.0)
# Paso de temperatura
dT = np.float(0.1)
# Agrego el detalle cerca de la temperatura crítica
detalle = [2.15,2.17,2.23, 2.25, 2.27, 2.32, 2.35, 2.45] 
# Número de pasos para la primer corrida (termalización)
N_term = '50000'
# Número de pasos para la segunda corrida (medición)
N_medi = '1000000'
# Número de corridas para cada temperatura (múltiplo del # de cores)
Nrun = 20
#
# FIN PARAMETROS DE ENTRADA
###############################################################################

# Lista de temperaturas
tempe = np.arange(T_min,T_max+dT,dT)
tempe = tempe.tolist() + detalle
tempe.sort(reverse=True)
# lo convierto a string
tempe = [str(i) for i in tempe]

# Cantidad de corridas por core
Nrun_local = Nrun//size
# Llista de las corridas para cada core
run_local = range(rank*Nrun_local,(rank+1)*Nrun_local)
aviso = np.ones(1)

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()

###############################################################################
# LOOP QUE CORRE SOBRE TODAS LAS TEMPERATURAS
###############################################################################
T_anterior = []    # Buffer para copiar el estado final a T anterior
                   # al estado actual.
for T in tempe:
    # Nombre de la carpeta uqe se va a crear
    carpeta = T + '_tmpfolder'
    # Camino completo de la carpeta que se va a crear
    path_carpeta = os.path.join(curr_dir,carpeta)
    # Sólo un core se encarga de armar el directorio de temperatura
    if rank==0:
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
        isf.escribe_entrada('T',T)
        isf.escribe_entrada('N',N_term)
        # Sólo para bloquear al resto de los procesos hasta que el root haya
        # hecho la carpeta. No sé cómo hacerlo más directo.
        for m in range(1,size):
            comm.Send(aviso,dest=m)
        print('Core {0} ya armó el directorio para T={1}'.format(rank,T))
    else:
        comm.Recv(aviso,source=0) # Bloquea al resto hasta recibir mensaje
        os.chdir(path_carpeta)    # Se meten en la carpeta ya creada por root
        
    ###########################################################################
    # DISTINTAS CORRIDAS A LA MISMA T PARA OBTENER ERROR ESTADISTICO
    #      - LOOP PARALELIZADO
    ###########################################################################
    for i in run_local:
    # Nombre de las carpetas con las corridas
        carpeta_runs = 'RUN' + '%02i'%i
        # Camino de las carpetas con las corridas
        path_runs = os.path.join(path_carpeta,carpeta_runs)
        # Crea las carpetas con las corridas
        os.makedirs(path_runs)
        # Copia archivos a cada carpeta
        shutil.copy('parametros.dat',path_runs)
        os.chdir(path_runs)
        # Escribe la semmilla en la carpeta
        isf.escribe_semilla()

        #########################################################
        ######### Para utilizar el estado de temperatura anterior
        isf.copia_estado_temp_anterior(path_runs,T_anterior,T)
        #########################################################
        
        #######################################################################
        # EJECUTA EL PROGRAMA ISING - TERMALIZACIÓN
        #######################################################################
        # Esto funciona para python >= 2.7         
        #salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/ising'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        #print(salida)        
        # Guardo la salida para ver ue hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        # Guardo las salidas por si hacen falta
        os.rename('energia.dat','energia_terma.dat')
        os.rename('magneti.dat','magneti_terma.dat')    

        #########################################################
        ######### Para utilizar el estado de temperatura anterior
        os.rename('estado.dat','ultimo_estado.dat')
        #########################################################
        
        #######################################################################
        # EJECUTA EL PROGRAMA ISING - MEDICIÓN
        #######################################################################
        # Corre por segunda vez tomando el estado anterior. Aumento el N
        isf.escribe_entrada('N',N_medi)
        print('Core {0} corriendo {1} a la temperatura {2}'.format(rank,carpeta_runs,T))
        # Esto funciona para python >= 2.7           
        # salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/ising'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]  
        # Guardo la salida para ver ue hizo
        with open('log2.txt','w') as arch: arch.write(salida) 
        # Copia los valores medios
        isf.copia_val_medios_runs(i,path_carpeta)
        # Sale de la carpeta de corridas        
        os.chdir(path_carpeta)
    # Espera a que todos los cores teminen sus procesos
    comm.Barrier()
    #
    # FIN DEL LOOP PARALELIZADO
    ###########################################################################
    if rank==0:
        # Hace estadística de todas las corridas y lo copia en un archivo
        isf.copia_val_medios(T,Nrun,curr_dir)
        # Sale de la carpeta
        os.chdir(curr_dir)        
    # Guardo la temperatura para copiar archivos luego
    T_anterior = T

###############################################################################
