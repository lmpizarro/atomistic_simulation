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
import dinmol_funcs as dm

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
if os.path.isfile("parametros_t.py"):
    
    import parametros_t as parametros
  
    N_part   = parametros.N_part
    Rho      = np.float(parametros.Rho)
    Temp_min = np.float(parametros.Temp_min)
    Temp_max = np.float(parametros.Temp_max)
    dTemp    = np.float(parametros.dTemp)
    #T_detail_min = np.float(parametros.T_detail_min)
    #T_detail_max = np.float(parametros.T_detail_max)
    #dT_detail = np.float(parametros.dT_detail)
    dt       = np.float(parametros.dt)
    N_term   = parametros.N_term
    N_medi   = parametros.N_medi
    N_grab   = parametros.N_grab
    Nrun     = parametros.Nrun
else:
    # Número de partículas
    N_part = 200
    # Densidad de partículas
    Rho = 0.3
    #------ Barrido de temperaturas
    # Temperatura mínima
    Temp_min = 0.7
    # Temperatura máxima
    Temp_max = 1.4
    # Paso de temperatura
    dTemp = 0.05
    # Agrego el detalle cerca de la temperatura crítica
    #T_detail_min = 2.10
    #T_detail_max = 2.50
    #dT_detail = 0.02
    # abs(K_grab) Cada cuántos puntos se quiere grabar el archivo temporal
    # K_grab < 0 especifica que no se grabe ningún archivo temporal
    N_grab = 10
    # Paso temporal de integración
    dt = 0.001
    # Número de pasos para la primer corrida (termalización)
    N_term = '4000'
    # Número de pasos para la segunda corrida (medición)
    N_medi = '10000'
    # Número de corridas para cada temperatura
    Nrun = 8

# Calcula el lado del cubo para la densidad y número de partículas especificado
L = (N_part/Rho)^(1/3)
# Escribe los valores al archivo parametros.dat
# Lo hace el root
if rank==0:
    # Escribe la cantidad de partículas
    dm.escribe_entrada('N_part',str(N_part))
    # Escribe el lado del cubo
    dm.escribe_entrada('L',str(L))
    # Escribe el dt de la integración temporal
    dm.escribe_entrada('dt',str(dt))
    # Escribe cada cuánto se hacen mediciones y si se guardan archivos
    dm.escribe_entrada('N_medi',str(N_medi))

# Todos esperan a que root haya terminado
comm.Barrier()
#
# FIN PARAMETROS DE ENTRADA
###############################################################################

# Lista de temperaturas de la zona de detalle
detalle = np.arange(T_detail_min,T_detail_max + dT_detail, dT_detail)

# Lista de temperaturas de paso grueso
tempe = np.arange(T_min,T_max+dT,dT)

tempe = np.append(tempe, detalle)
# filtra por valores unicos de temperatura
# se debe redondear (con algún criterio) para evitar falsas duplicaciones
# por errores de redondeo de python
tempe = np.unique(tempe.round(decimals=6))
# ordena de mayor a menor el array numpy
tempe = np.fliplr([tempe])[0]

# Cantidad de corridas por core 
Nrun_local = Nrun//size

# Llista de las corridas para cada core
run_local = range(rank*Nrun_local,(rank+1)*Nrun_local)

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()

###############################################################################
# LOOP QUE CORRE SOBRE TODAS LAS TEMPERATURAS
###############################################################################
T_anterior = []    # Buffer para copiar el estado final a T anterior
                   # al estado actual.
for T in tempe:
    Tnombre = str(T)
    # Nombre de la carpeta uqe se va a crear
    carpeta = Tnombre + '_tmpfolder'
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
                print('Ya existe el directorio. Se omite corrida con T=' + Tnombre)
                T_anterior = Tnombre
                # Si se llega al final de la lista, se finaliza el programa
                if T==tempe[-1]:
                    print('No hay nada por hacer, se sale del programa')
                    comm.Abort()
                continue
        # Copia el archivo de entrada a la carpeta 
        shutil.copy('parametros.dat',path_carpeta)
        # Se mete en la carpeta
        os.chdir(path_carpeta)       
        # Cambia el archivo de entrada adentro de la carpeta
        isf.escribe_entrada('T',Tnombre)
        isf.escribe_entrada('N',N_term)
        print('Core {0} ya armó el directorio para T={1}'.format(rank,T))
    # Una vez que root encontró la carpeta donde debe trabajar, manda la información
    # al resto de los procesos
    T            = comm.bcast(T, root=0)
    T_anterior   = comm.bcast(T_anterior, root=0)
    Tnombre      = comm.bcast(Tnombre, root=0)
    path_carpeta = comm.bcast(path_carpeta, root=0)
    # Se meten en la carpeta
    os.chdir(path_carpeta) 
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
        isf.copia_estado_temp_anterior(path_runs,T_anterior,Tnombre)
        #########################################################
        
        #######################################################################
        # EJECUTA EL PROGRAMA ISING - TERMALIZACIÓN
        #######################################################################
        # Esto funciona para python >= 2.7         
        #salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/ising'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        # Guardo las salidas por si hacen falta
        if os.path.isfile('energia.dat'):
            os.rename('energia.dat','energia_terma.dat')
        elif os.path.isfile('magneti.dat'):
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
        # Lee datos del archivo <val_medios.dat> y <aceptaciones.dat>
        datos = isf.lee_datos_runs(i)
        # Root se encarga de recibir los datos y procesarlos
        if rank==0:
            # Root escribe sus propios datos en <runs_esatadistica.dat>
            isf.escribe_datos_runs(datos,path_carpeta) 
            for j in range(1,size):
                # root recibe los datos del resto de los procesos
                datos_recv = comm.recv(source=MPI.ANY_SOURCE)
                # root escribe los datos en <runs_estadistica.dat>
                isf.escribe_datos_runs(datos_recv,path_carpeta) 
        else:
            # El resto de los procesos le mandan los datos a root
            comm.send(datos,dest=0)
        # Sale de la carpeta de corridas        
        os.chdir(path_carpeta)
    #
    # FIN DEL LOOP PARALELIZADO
    ###########################################################################
    if rank==0:
        # Hace estadística de todas las corridas y copia en <tablas_temperatura.dat> 
        isf.copia_val_medios(T,Nrun,curr_dir)
        # Sale de la carpeta
        os.chdir(curr_dir)        
    # Guardo la temperatura para copiar archivos luego
    T_anterior = Tnombre

###############################################################################
