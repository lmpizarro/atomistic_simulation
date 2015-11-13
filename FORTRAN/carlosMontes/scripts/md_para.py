#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import time

import numpy as np 
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#/home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -np 4 /home/pizarro/python/bin/python ./Melting.py

# inicializa la semilla del generador de números aleatorios
# de python para crear seed.dat 
np.random.seed(int(time.time()))
# para tener el máximo y el mínimo entero
ii32 = np.iinfo(np.int32)

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
    print "El archivo param_md_para.py debe estar inicializado. "
    print "El programa terminará inmediatamente "
    exit()

# Calcula el lado del cubo para la densidad y número de partículas especificado
L = np.power( N_part/Rho , np.float(1)/3 ) 


# Directorio raíz donde está el ejecutable y este script
root_dir = os.getcwd()

# Lista de temperaturas de paso grueso
Temperatures = np.arange(Temp_min,Temp_max+dTemp,dTemp)
# Ordeno de mayor a menor las temperaturas
Temperatures = np.fliplr([Temperatures])[0]

Cant_Puntos = Temperatures.size
Puntos_Por_Nodo = (Cant_Puntos / size) + 1

Ti = rank * Puntos_Por_Nodo
Tf = (rank + 1)  * Puntos_Por_Nodo
temps =  Temperatures[Ti: Tf]

print "el nodo: ", rank, "calcula las temps", temps
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

def crea_seed_dat ():
    semilla = np.random.randint(ii32.min, ii32.max)
    with open("seed.dat",'w') as arch: arch.write(str(semilla) + "\n")
#
# Corrida cero (0)
#
if rank == 0:
    corrida_cero_path = "corrida0" 
    if not os.path.exists(corrida_cero_path):
        os.makedirs(corrida_cero_path)
    os.chdir(corrida_cero_path)       
    gen_parameters_dat(Temperatures[0], N_term)
    
    crea_seed_dat()

    proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
    salida = proc.communicate()[0]
    # Guardo la salida para ver que hizo
    with open('log1.txt','w') as arch: arch.write(salida)
    with open('corrida0.txt','w') as arch: arch.write("corrida0" + " " + str(Temperatures[0]))
    os.chdir("../")
# espera a que el nodo 0 termine de calcular esta DM
comm.Barrier()

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

            crea_seed_dat()

            # Copia el archivo de estados  de corrida0 a la carpeta de temperatura
            shutil.copy('../../corrida0/estados.dat',"./")
            os.chdir("../")
        os.chdir("../")
# espera a que el nodo 0 termine de crear las carpetas necesarias
comm.Barrier()

#
# Corre la termalización
# y la post termalización
for iteration in range(Nrun):
    iteration_path = "iteration"+ "_" + str(iteration)
    os.chdir(iteration_path)

    for ts in temps:
        carpeta_temp = "temp" + "_" + str(ts)
        os.chdir(carpeta_temp)
        # corrida de la DM 
        proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log1.txt','w') as arch: arch.write(salida)
        with open('termalizacion.txt','w') as arch: arch.write("termalizado" + " " + str(ts))
        # Se crean los parámetros para las corridas 
        gen_parameters_dat(ts, N_medi)

        proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        # Guardo la salida para ver que hizo
        with open('log2.txt','w') as arch: arch.write(salida)

        os.chdir("../")       
    os.chdir("../")       
#
# Corre luego de la termalización 
#
#for iteration in range(Nrun):
#    iteration_path = "iteration"+ "_" + str(iteration)
#    os.chdir(iteration_path)
#
#    for ts in temps:
#        carpeta_temp = "temp" + "_" + str(ts)
#        os.chdir(carpeta_temp)
#
#        proc = subprocess.Popen([root_dir+'/dinmod'],stdout=subprocess.PIPE)
#        salida = proc.communicate()[0]
#        # Guardo la salida para ver que hizo
#        with open('log2.txt','w') as arch: arch.write(salida)
#
#        os.chdir("../")       
#    os.chdir("../")

if rank == 0:
    print "fin: ", rank 
