# -*- coding: utf-8 -*-

#!/usr/bin/python

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
###############################################################################

import os
import os.path
import shutil
import subprocess
import errno
import numpy as np

# Funciones propias para este script
import ising_fun as isf

###############################################################################       
#   PARAMETROS DE ENTRADA
###############################################################################

if os.path.isfile("parametros.py"):
   
    import parametros
    
    N_red = parametros.N_red
    M_red = parametros.M_red
    T_min = np.float(parametros.T_min)
    T_max = np.float(parametros.T_max)
    dT = np.float(parametros.dT)
    T_detail_min = np.float(parametros.T_detail_min)
    T_detail_max = np.float(parametros.T_detail_max)
    dT_detail = np.float(parametros.dT_detail)
    N_term  = parametros.N_term
    N_medi  = parametros.N_medi
    Nrun = parametros.Nrun
    N_grab = parametros.N_grab
else:
    ## Tamaño de la red de spines
    N_red = 20
    M_red = 20
    # Cada cuántos puntos se quiere grabar el archivo temporal
    N_grab = 0 
   
    # Barrido de temperaturas
    # Temperatura mínima
    T_min = np.float(1.5)
    # Temperatura máxima
    T_max = np.float(1.7)
    # Paso de temperatura
    dT = np.float(0.1)
    # Agrego el detalle cerca de la temperatura crítica
    T_detail_min = np.float(2.15)
    T_detail_max = np.float(2.45)
    dT_detail = np.float(0.03)
    # Número de pasos para la primer corrida (termalización)
    N_term = '4000'
    # Número de pasos para la segunda corrida (medición)
    N_medi = '10000'
    # Número de corridas para cada temperatura
    Nrun = 10
    
# Escribe los valores al archivo parametros.dat
# Especifica tamaño de la red de spines
isf.escribe_entrada('N_red',str(N_red))
isf.escribe_entrada('M_red',str(M_red))
# Especifica cada cuanto se graban datos temporales 
isf.escribe_entrada('N_grab',str(N_grab))

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

# Número de corridas para las T fuera de la zona crítica
N_run_base = Nrun

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
    # Se crea la carpeta sólo si ésta no existe. De lo contrario se saltea la 
    # temperatura correspondiente
    try: 
        os.makedirs(path_carpeta)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print('Ya existe el directorio. Se omite corrida con T=' + Tnombre)
            continue
    # Copia el archivo de entrada a la carpeta 
    shutil.copy('parametros.dat',path_carpeta)
    #shutil.copy('ising',path_carpeta) 
    # Se mete en la carpeta
    os.chdir(path_carpeta)       
    # Cambia el archivo de entrada adentro de la carpeta
    isf.escribe_entrada('T',Tnombre)
    # Corre el programa para ver la convergencia
    isf.escribe_entrada('N',N_term) 

    ##########################################################################
    # DISTINTAS CORRIDAS A LA MISMA T PARA OBTENER ERROR ESTADISTICO
    ##########################################################################

    # Determino si estoy en la zona crítica
    # para aumentar la cantidad de corridas
    if T >= T_detail_min and T <= T_detail_max:
        Nrun = 2 * N_run_base
        print "Zona crítica: ", Nrun
    else:
        Nrun = N_run_base
        print "Zona No crítica: ", Nrun

    for i in range(0,Nrun):
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
        # salida = subprocess.check_output(curr_dir+'/ising')
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
        print('Corriendo {0} a la temperatura {1}'.format(carpeta_runs,T))
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
    # Hace estadística de todas las corridas y lo copia en un archivo
    isf.copia_val_medios(T,Nrun,curr_dir)
    # Sale de la carpeta
    os.chdir(curr_dir)
    T_anterior = Tnombre
    
##############################################################################

