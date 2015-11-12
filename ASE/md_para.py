#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import time

import numpy as np 
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

from asap3 import Atoms, EMT, units
from ase.visualize.primiplotter import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.langevin import Langevin

from ase.io.trajectory import Trajectory

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# /home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -n 4 python md_para.py

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
    N_run     = np.int(parametros.Nrun)
    epsilon = np.float(parametros.epsilon)
    sigma   = np.float(parametros.sigma)
    masa = np.float(parametros.masa)
    N_save  = np.int(parametros.N_save)
    gamma   = np.float(parametros.gamma)
else:
    print "El archivo param_md_para.py debe estar inicializado. "
    print "El programa terminará inmediatamente "
    exit()

# Directorio raíz donde está el ejecutable y este script
root_dir = os.getcwd()

# Lista de temperaturas de paso grueso
Temperatures = np.arange(Temp_min, Temp_max + dTemp, dTemp)
# Ordeno de mayor a menor las temperaturas
Temperatures = np.fliplr([Temperatures])[0]

Cant_Puntos = Temperatures.size
Puntos_Por_Nodo = (float(Cant_Puntos) / float(size)) 

Ti = rank * Puntos_Por_Nodo
Tf = (rank + 1)  * Puntos_Por_Nodo
temps =  Temperatures[Ti: Tf]

print rank,  Ti, Tf

if rank == 0:
    print "nro de iteraciones ", N_run 
    print "Temperatures: ", Temperatures
    print Temperatures[0], Temperatures[Temperatures.size -1], Temperatures.size
comm.Barrier()

print "el nodo: ", rank, "calcula las temps", temps

# Create the atoms
atoms = FaceCenteredCubic(size = parametros.Size, 
                          symbol="Cu", 
                          pbc=parametros.Pbc,
                          latticeconstant = parametros.LatticeConstant)

# Associate the EMT potential with the atoms
atoms.set_calculator(EMT())


# Make the Langevin dynamics module
dyn = Langevin(atoms, 5*units.fs, units.kB*Temperatures[0], 0.002)

# Some functions for calculating the actual temperature, energy, ...

def temperature(a):
    return 2.0/3.0 * a.get_kinetic_energy() / (len(a) * units.kB)

def etotal(a):
    return (a.get_kinetic_energy() + a.get_potential_energy()) / len(a)

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

    with open('corrida0.txt','w') as arch: arch.write("corrida0" + " " + str(Temperatures[0]))
    os.chdir("../")
# espera a que el nodo 0 termine de calcular esta DM
comm.Barrier()

#
# Creación de las carpetas
# y parametros.dat en cada directorio
# para termalizar
for iteration in range(N_run):
    if rank == 0:
        iteration_path = "iteration"+ "_" + str(iteration + 1)
        if not os.path.exists(iteration_path):
            os.makedirs(iteration_path)
        os.chdir(iteration_path)       
        for ts in Temperatures:
            carpeta_temp = 'temp_%04d' % (ts)
            if not os.path.exists(carpeta_temp):
               os.makedirs(carpeta_temp)

        os.chdir("../")
# espera a que el nodo 0 termine de crear las carpetas necesarias
comm.Barrier()


#
# Corre la termalización
# y la post termalización
for iteration in range(N_run):
    iteration_path = "iteration"+ "_" + str(iteration + 1)
    os.chdir(iteration_path)

    for ts in temps:
        carpeta_temp = 'temp_%04d' % (ts)
        os.chdir(carpeta_temp)
        # corrida de la DM 
        dyn.set_temperature(units.kB*ts)

        with open('calcs.txt','w') as arch: arch.write("Calculos a temperatura: " + str(ts) + "\n")
        for i in range(N_medi/100):
            dyn.run(100)
            atoms_file_name = 'atoms_%04d.xyz' % (i)
            atoms.write(atoms_file_name)
            message_calc = "%-10.5f  %.0f \n" % (etotal(atoms), temperature(atoms))
            with open('calcs.txt','a+') as arch: arch.write(message_calc)
            print (("%04d %04d %.0f %04d %s")%(iteration, rank, ts, i, atoms_file_name))

        mess_log = "corrida iteracion: " + str(iteration + 1) + " temperatura: " + str(ts)
        with open('log2.txt','w') as arch: arch.write(mess_log)

        os.chdir("../")       
    os.chdir("../")
 
if rank == 0:
    print "fin: ", rank 
