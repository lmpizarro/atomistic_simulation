#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

import numpy as np 
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#/home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -np 4 /home/pizarro/python/bin/python ./Melting.py


###############################################################################       
#   PARAMETROS DE ENTRADA
###############################################################################
if os.path.isfile("param_md_para.py"):
    
    import param_md_para as parametros
  
    N_part   = np.int(parametros.N_part)
    Rho      = np.float(parametros.Rho)
    Temp_min = np.float(parametros.Temp_min)
    Temp_max = np.float(parametros.Temp_max)
    dTemp    = np.float(parametros.dTemp)
    #T_detail_min = np.float(parametros.T_detail_min)
    #T_detail_max = np.float(parametros.T_detail_max)
    #dT_detail = np.float(parametros.dT_detail)
    dt       = np.float(parametros.dt)
    N_term   = np.int(parametros.N_term)
    N_medi   = np.int(parametros.N_medi)
    N_grab   = np.int(parametros.N_grab)
    Nrun     = np.int(parametros.Nrun)
    epsilon = np.float(parametros.epsilon)
    sigma   = np.float(parametros.sigma)
    masa = np.float(parametros.masa)
    N_save  = np.int(parametros.N_save)
    gamma   = np.float(parametros.gamma)

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
    N_term = 4000
    # Número de pasos para la segunda corrida (medición)
    N_medi = 10000
    # Número de corridas para cada temperatura
    Nrun = 8
    epsilon = 1.0
    sigma   = 1.0
    masa = 1.0
    N_save  = 1000
    gamma   = 1.0

# Directorio raíz donde está el ejecutable y este script
root_dir = os.getcwd()

# Calcula el lado del cubo para la densidad y número de partículas especificado
L = np.power( N_part/Rho , np.float(1)/3 ) 

# Lista de temperaturas de paso grueso
Temperatures = np.arange(Temp_min,Temp_max+dTemp,dTemp)

Cant_Puntos = Temperatures .size
Puntos_Por_Nodo = Cant_Puntos / size + 1

Ti = rank * Puntos_Por_Nodo
Tf = (rank + 1)  * Puntos_Por_Nodo
temps =  Temperatures[Ti: Tf]

print rank, temps
#
# Formato de parametros.dat
# Temp Npart L tinteg pasadas sigma epsilon masa Ngrabacion
# Gamma
#1.0 200 5.5032 0.001 20000 1.0 1.0 1.0 100
#0.5
#
def gen_parameters_dat (temp, Nmedi):
    file_name = "parametros.dat"
    p_corr = str(temp) + " " + str(N_part) + " " + str(L) + " " + str(dt) + " " + str(Nmedi)   
    p_part = str(sigma) + " " + str(epsilon) + " " + str(masa)
    str_to_write = p_corr + " " + p_part +  " " + str(N_save) + "\n" + str(gamma)
    with open(file_name,'w') as arch: arch.write(str_to_write)
    #fo = open(name, "wb")
    #fo.write(str_to_write)
    #fo.write(str(sigma))
    #fo.close()

def gen_data (rank):
    name = "foo"+ str(rank)  + ".txt"
    fo = open(name, "wb")
    fo.write("hola mundo: " + str(rank) )
    fo.close()

#
# Creación de las carpetas
# y parametros.dat en cada directorio
# para termalizar
for iteration in range(Nrun):
    if rank == 0:
        iteration_path = "iteration"+ "_" + str(iteration)
        if not os.path.exists(iteration_path):
            os.makedirs(iteration_path)
        os.chdir(iteration_path)       
        for ts in Temperatures:
            carpeta_temp = "temp" + "_" + str(ts)
            if not os.path.exists(carpeta_temp):
               os.makedirs(carpeta_temp)
            os.chdir(carpeta_temp)
            # Se crean los parámetros para la termalización
            gen_parameters_dat(ts, N_term)
            os.chdir("../")
        os.chdir("../")
comm.Barrier()

#
# Corre la termalización
#
for iteration in range(Nrun):
    iteration_path = "iteration"+ "_" + str(iteration)
    os.chdir(iteration_path)

    for ts in temps:
        carpeta_temp = "temp" + "_" + str(ts)
        os.chdir(carpeta_temp)
        proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        with open('termalizacion.txt','w') as arch: arch.write("termalizado" + " " + str(ts))
        # Se crean los parámetros para las corridas 
        gen_parameters_dat(ts, N_medi)
        os.chdir("../")       

    os.chdir("../")       


#
# Corre luego de la termalización 
#
for iteration in range(Nrun):
    iteration_path = "iteration"+ "_" + str(iteration)
    os.chdir(iteration_path)

    for ts in temps:
        carpeta_temp = "temp" + "_" + str(ts)
        os.chdir(carpeta_temp)
        proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        # Se crean los parámetros para las corridas 
        os.chdir("../")       

    os.chdir("../")       






if rank == 0:
    print "fin: ", rank 

