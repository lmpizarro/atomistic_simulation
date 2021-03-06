#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################       
#  SE EJECUTA CON: mpirun.mpich -n 4  python corridas_paralelo_rho.py
###############################################################################
#
# Script para correr el programa de ising a distinta densidades
# Crea una carpeta para cada densidad. En cada una de esas carpetas, a su
# vez, crea Nruns carpetas para hacer estadística y obtener los valores con
# sus respectivos errores.
#
# En cada densidad, se utiliza el valor final de la densidad anterior.
# Arbitrariamente, se toma el valor final de RUN00 como el valor inicial de
# todas las corridas a la siguiente temperatura (menor)
#
# Los resultados los guarda en el archivo 'tabla_densidad.dat'
# Cada fila de 'tabla_densidad.dat' representa los resultados a una dada
# densidad con sus respectivos errores. Las columnas son:
# T rho <p> std(p) <sp> std(sp) <Temp> std(Temp) <U> std(U)
# T es la temperatura fijada. Temp es la medida.
# Los errores se calculan como std(RUN)/sqrt(NRUN)
#
# Se paralelizan las corridas que se hace a una dada densidad. La limitación
# está en que se debe esperar que terminen las corridas a la densidad anterior
# para utilizar el estado final. De lo contrario se podría parelelizar también
# la corrida para cada densidad (en el archivo 'sin errores' está hecho así)
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
if os.path.isfile("parametros_rho.py"):
    
    import parametros_rho as parametros
  
    N_part   = np.int(parametros.N_part)
    Temp     = np.float(parametros.Temp)
    Rho      = np.asarray(parametros.Rho)
    #Rho_min  = np.float(parametros.Rho_min)
    #Rho_max  = np.float(parametros.Rho_max)
    #dRho     = np.float(parametros.dRho)
    Rho_detail_min = np.float(parametros.Rho_detail_min)
    Rho_detail_max = np.float(parametros.Rho_detail_max)
    dRho_detail = np.float(parametros.dRho_detail)
    dt       = np.float(parametros.dt)
    N_term   = np.int(parametros.N_term)
    N_medi   = np.int(parametros.N_medi)
    N_grab   = np.int(parametros.N_grab)
    Nrun     = np.int(parametros.Nrun)
else:
    # Número de partículas
    N_part = 200
    # Temperatura 
    Temp = 1.1
    #------ Barrido de densidades
    # Densidad mínima
    #Rho_min = 0.7
    # Densidad máxima
    # Rho_max = 1.4
    # Paso de densidad 
    # dRho = 0.05
    # Agrego el detalle cerca de la temperatura crítica
    Rho_detail_min = 2.10
    Rho_detail_max = 2.50
    dRho_detail = 0.02
    # abs(K_grab) Cada cuántos puntos se quiere grabar el archivo temporal
    # K_grab < 0 especifica que no se grabe ningún archivo temporal
    N_grab = 10
    # Paso temporal de integración
    dt = 0.001
    # Número de pasos para la primer corrida (termalización)
    N_term = 4000
    # Número de pasos para la segunda corrida (medición)
    N_medi = 10000
    # Número de corridas para cada temperatura
    Nrun = 8

# Escribe los valores al archivo parametros.dat
# Lo hace el root
if rank==0:
    # Escribe la cantidad de partículas
    dm.escribe_entrada('N_part',str(N_part))
    # Escribe la temperatura 
    dm.escribe_entrada('Temp',str(Temp))

# Todos esperan a que root haya terminado
comm.Barrier()
#
# FIN PARAMETROS DE ENTRADA
###############################################################################

# Lista de temperaturas de la zona de detalle
detalle = np.arange(Rho_detail_min,Rho_detail_max + dRho_detail, dRho_detail)

# Lista de densidades de paso grueso
Rho = np.array(Rho)

Rho = np.append(Rho, detalle)
# filtra por valores unicos de densidad 
# se debe redondear (con algún criterio) para evitar falsas duplicaciones
# por errores de redondeo de python
Rho = np.unique(Rho.round(decimals=6))
# ordena de mayor a menor el array numpy
Rho = np.fliplr([Rho])[0]

# Cantidad de corridas por core
# Notar que Nrun conviene que sea múltiplo de la cantidad
# de procesos con que se corra, de lo contrario se redondea para abajo.
Nrun_local = Nrun//size

# Llista de las corridas para cada core
run_local = range(rank*Nrun_local,(rank+1)*Nrun_local)

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()

###############################################################################
# CORRIDA INICIAL - MINIMIZACIÓN DE LA ENERGÍA 
###############################################################################
# Primera corrida para separar las partículas minimizando la energía potencial.
# Lo corre sólo root y luego se copia el archivo al resto de las carpetas

if rank==0:
    # Escribe la densidad inicial (la mayor), el L más chico
    dm.escribe_entrada( 'L',dm.rho_2_lado( N_part, Rho[0] ) )
    # Escribe el número de pasos temporale
    dm.escribe_entrada('N_pasos','3000')
    # Escribe el dt de la integración temporal
    dm.escribe_entrada('dt','0.0001')
    # Escribe cada cuánto se hacen mediciones y si se guardan archivos
    dm.escribe_entrada('N_grab','1')
        
    print('Core {0} corriendo la minimización de energía a  Rho={1}'.format(rank,Rho[0]))
    
    proc = subprocess.Popen(['./dinmod'],stdout=subprocess.PIPE)
    salida = proc.communicate()[0]
    # Guardo la salida para ver que hizo
    with open('log0.txt','w') as arch: arch.write(salida)
    
    # Reescribe parametros para dejar listo el archivo de entrada
    # Escribe número de pasos para la termalización
    dm.escribe_entrada('N_pasos',str(N_term))
    # Escribe dt
    dm.escribe_entrada('dt',str(dt))
    # Escribe cada cuánto se hacen mediciones y si se guardan archivos
    dm.escribe_entrada('N_grab',str(N_grab))

# Todos esperan a que root haya terminado
comm.Barrier()

###############################################################################
# LOOP QUE CORRE SOBRE TODAS LAS DENSIDADES 
###############################################################################
R_anterior = []  # Buffer para copiar el estado final a R anterior
                   # al estado actual.
for R in Rho:
    Rnombre = str(R)
    # Nombre de la carpeta uqe se va a crear
    carpeta = 'densidad_' + Rnombre 
    # Camino completo de la carpeta que se va a crear
    path_carpeta = os.path.join(curr_dir,carpeta)
    # Sólo un core se encarga de armar el directorio de densidad
    if rank==0:
        # Se crea la carpeta sólo si ésta no existe. De lo contrario se saltea la 
        # densidad correspondiente
        try: 
            os.makedirs(path_carpeta)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            else:
                print('Ya existe el directorio. Se omite corrida con Rho=' + Rnombre)
                R_anterior = Rnombre
                # Si se llega al final de la lista, se finaliza el programa
                if R==Rho[-1]:
                    print('No hay nada por hacer, se sale del programa')
                    comm.Abort()
                continue
        # Copia el archivo de estados la carpeta de densidad 
        shutil.copy('estados.dat',path_carpeta)
        # Copia el archivo de entrada a la carpeta de densidad
        shutil.copy('parametros.dat',path_carpeta)
        # Se mete en la carpeta
        os.chdir(path_carpeta)
        # Cambia el archivo de entrada adentro de la carpeta
        dm.escribe_entrada('L', dm.rho_2_lado(N_part,R))
        # para utilizar el estado de temperatura anterior
        dm.copia_estado_temp_anterior(path_carpeta,R_anterior,Rnombre)
        print('Core {0} ya armó el directorio para Rho={1}'.format(rank,R))
    # Una vez que root encontró la carpeta donde debe trabajar, manda la información
    # al resto de los procesos
    R            = comm.bcast(R, root=0)
    R_anterior   = comm.bcast(R_anterior, root=0)
    Rnombre      = comm.bcast(Rnombre, root=0)
    path_carpeta = comm.bcast(path_carpeta, root=0)
    # Se meten en la carpeta
    os.chdir(path_carpeta) 
    
    ###########################################################################
    # DISTINTAS CORRIDAS A LA MISMA Rho PARA OBTENER ERROR ESTADISTICO
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
        # Copia archivos de estados a cada carpeta
        shutil.copy('estados.dat',path_runs)
        # Pasa a carmeta de cada corrida RUN..
        os.chdir(path_runs)
        # Escribe la semmilla en la carpeta
        dm.escribe_semilla()

        #########################################################
        ######### Para utilizar el estado de temperatura anterior
        #dm.copia_estado_temp_anterior(path_runs,T_anterior,Tnombre)
        #########################################################
        
        #######################################################################
        # EJECUTA EL PROGRAMA DINMOD - TERMALIZACIÓN
        #######################################################################
        # Esto funciona para python >= 2.7         
        #salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        # Guardo las salidas por si hacen falta
        if os.path.isfile('energias.dat'):
            os.rename('energias.dat','energias_terma.dat')
        elif os.path.isfile('presion.dat'):
            os.rename('presion.dat','presion_terma.dat')    
        elif os.path.isfile('temperatura.dat'):
            os.rename('temperatura.dat','temperatura_terma.dat')
        
        #######################################################################
        # EJECUTA EL PROGRAMA ISING - MEDICIÓN
        #######################################################################
        # Corre por segunda vez tomando el estado anterior. Aumento el N
        dm.escribe_entrada('N_pasos',str(N_medi))
        print('Core {0} corriendo {1} a la densidad {2}'.format(rank,carpeta_runs,R))
        # Esto funciona para python >= 2.7           
        # salida = subprocess.check_output(curr_dir+'/dinmod')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]  
        # Guardo la salida para ver ue hizo
        with open('log2.txt','w') as arch: arch.write(salida) 
        # Lee datos del archivo <val_medios.dat> 
        datos,header = dm.lee_datos_runs(i)
        # Root se encarga de recibir los datos y procesarlos
        if rank==0:
            # Root escribe sus propios datos en <runs_esatadistica.dat>
            dm.escribe_datos_runs(datos,header,path_carpeta) 
            for j in range(1,size):
                # root recibe los datos del resto de los procesos
                datos_recv = comm.recv(source=MPI.ANY_SOURCE)
                # root escribe los datos en <runs_estadistica.dat>
                dm.escribe_datos_runs(datos_recv,header,path_carpeta) 
        else:
            # El resto de los procesos le mandan los datos a root
            comm.send(datos,dest=0)
        # Sale de la carpeta de corridas        
        os.chdir(path_carpeta)
    #
    # FIN DEL LOOP PARALELIZADO
    ###########################################################################
    if rank==0:
        # Hace estadística de todas las corridas y copia en <tabla_presion.dat> 
        dm.copia_val_medios(Temp,R,Nrun,curr_dir,'tabla_presion.dat')
        # Sale de la carpeta
        os.chdir(curr_dir)        
    # Guardo la temperatura para copiar archivos luego
    R_anterior = Rnombre

###############################################################################
