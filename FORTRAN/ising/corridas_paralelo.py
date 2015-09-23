#!/usr/bin/python
# -*- coding: utf-8 -*-
 
# Correr con:  mpirun.mpich -n 4  python corridas_paralelo.py 
# Crea una carpeta para cada temperatura. En cada una de esas carpetas, a su
# vez, crea Nruns carpetas para hacer estadística y obtener los valores con
# sus respectivos errores.

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

# Información de dónde está            
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
            exitshutil.copy
    with open('parametros.dat','w') as fileout:
        fileout.write(new)

def escribe_semilla():
    random.seed()
    # Obtengo un número aleatorio entero de 4 bytes
    entero_aleatorio = random.randint(-2147483648,2147483647)
    with open('seed.dat','w') as fileout:
        fileout.write(str(entero_aleatorio))
        
def copia_val_medios(T,N):
    E=[];M=[];c=[];s=[];a=[]
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('runs_estadistica.dat','r') as valruns:
        for line in valruns:
            col = [float(x) for x in line.split()]
            E.append(col[2])
            c.append(col[3])            
            M.append(col[4])
            s.append(col[5])
            a.append(col[6])
    arch_comun = os.path.join(curr_dir,'tablas_temperatura.dat')
    fa = 1/np.sqrt(N)
    print('Se promediaron {0} de {1} valores para T = {2}'.format(len(E),N,T))
    valor = [float(T), np.mean(E), np.std(E)*fa, np.mean(M), np.std(M)*fa,
                np.mean(c), np.std(c)*fa, np.mean(s), np.std(s)*fa,
                np.mean(a), np.std(a)*fa]
    # Abro el archivo donde se volcarán las resultados de cada temperatura   
    with open(arch_comun,'a') as comun:
        #valor_str = '  '.join( ('{:.5E}').format(j) for j in valor) 
        valor_str = '  '.join( format(j,'.5E') for j in valor) 
        comun.write( valor_str +'\n')   
            
def copia_val_medios_runs(i):
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('val_medios.dat','r') as vmed:
        valor=vmed.readline()
    # Lee el archivo con la cantidad de aceptaciones
    with open('aceptaciones.dat','r') as archacept:
        acep=archacept.readline()
        colu = [float(x) for x in acep.split()]
        # Guardará el porcentaje de aceptaciones
        por_acep = 100.0*colu[1]/colu[2]
    arch_comun = os.path.join(path_carpeta,'runs_estadistica.dat')
    with open(arch_comun,'a') as comun:
        comun.write(str(i)+' ' + valor.rstrip() + ' ' + str(por_acep) + '\n')
        
def copia_estado_temp_anterior(path_act,T_ant,T_act):
    # Copia y reescribe el último estado calculado a una temperatura anterior
    if T_ant != []:
        # Reemplazo sólo la parte de la temperatura en toda la ruta
        path_ant = path_act.replace(T_act,T_ant)
        # Ruta completa del archivo a copiar
        # arch_ant = os.path.join(path_ant,'estado.dat')
        # Para todas las corridas utiliza el ultimo estado de RUN0
        arch_ant = os.path.join(os.path.dirname(path_ant),'RUN0','estado.dat')        
        # Copia el archivo
        shutil.copy(arch_ant,path_act)
        # reescribe el archivo
        os.rename('estado.dat','ultimo_estado.dat')

##############################################################################
##############################################################################       
# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()
# Lista con las temperaturas que se desean calcular
T_min = np.float(0.5)
T_max = np.float(4.0)
dT    = np.float(0.1)
tempe = np.arange(T_min,T_max+dT,dT)
tempe = tempe.tolist() + [2.15,2.17,2.23, 2.25, 2.27, 2.32, 2.35, 2.45] # Mäs detalle en la Tc
tempe.sort(reverse=True)
# lo convierto a string
tempe = [str(i) for i in tempe]
# Número de pasos para la primer corrida (termalización)
N_term = '40000'
# Número de pasos para la segunda corrida (medición)
N_medi = '1000000'
# Número de corridas para cada temperatura
Nrun  = 12
# Cantidad de corridas por core
Nrun_local = Nrun//size
# Llista de las corridas para cada core
run_local = range(rank*Nrun_local,(rank+1)*Nrun_local)
aviso = np.ones(1)

T_anterior = []
# Loop para crear todos los directorios y correr el ejecutable en ellos
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
        escribe_entrada('T',T)
        escribe_entrada('N',N_term)
        # Sólo para bloquear al resto de los procesos hasta que el root haya
        # hecho la carpeta. No sé cómo hacerlo más directo.
        for m in range(1,size):
            comm.Send(aviso,dest=m)
        print('Core {0} ya armó el directorio para T={1}'.format(rank,T))
    else:
        comm.Recv(aviso,source=0) # Bloquea al resto hasta recibir mensaje
        os.chdir(path_carpeta) # Se meten en la carpeta ya creada por root
        
   # Loop para correr N veces con los mismos parámetros para calcular el error
    for i in run_local:
    # Nombre de las carpetas con las corridas
        carpeta_runs = 'RUN' + str(i)
        # Camino de las carpetas con las corridas
        path_runs = os.path.join(path_carpeta,carpeta_runs)
        # Crea las carpetas con las corridas
        os.makedirs(path_runs)
        # Copia archivos a cada carpeta
        shutil.copy('parametros.dat',path_runs)
        os.chdir(path_runs)
        # Escribe la semmilla en la carpeta
        escribe_semilla()
        ####################################
        ### Para utilizar el estado de temperatura anterior
        # Copio el último estado calculado al a temperatura anterior
        
        copia_estado_temp_anterior(path_runs,T_anterior,T)
        
        ####################################
        ####### Ejecuta el programa
        # Esto funciona para python >= 2.7         
        #salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/ising'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        #print(salida)        
        # Guardo la salida para ver ue hizo
        f=open("log1.txt",'w')
        f.write(salida)
        f.close() 
        # Guardo las salidas por si hacen falta
        os.rename('energia.dat','energia_terma.dat')
        os.rename('magneti.dat','magneti_terma.dat')    
        ####################################
        ### Para utilizar el estado de temperatura anterior
        # Acá tengo que reescribir el 'ultimo_estado' con lo que obtengo (estado)
        
        os.rename('estado.dat','ultimo_estado.dat')

        ####################################
        # Corre por segunda vez tomando el estado anterior. Aumento el N
        escribe_entrada('N',N_medi)
        print('Core {0} corriendo {1} a la temperatura {2}'.format(rank,carpeta_runs,T))
        # Esto funciona para python >= 2.7           
        # salida = subprocess.check_output(curr_dir+'/ising')
        # Alternativa para python 2.6        
        proc = subprocess.Popen([curr_dir+'/ising'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
    
        # Guardo la salida para ver ue hizo
        f=open("log2.txt",'w')
        f.write(salida)
        f.close()    
        # Copia los valores medios
        copia_val_medios_runs(i)
        # Sale de la carpeta de corridas        
        os.chdir(path_carpeta)
    # Espera a que todos los cores teminen sus procesos
    comm.Barrier()
    if rank==0:
        # Hace estadística de todas las corridas y lo copia en un archivo
        copia_val_medios(T,Nrun)
        # Sale de la carpeta
        os.chdir(curr_dir)        
    # Guardo la temperatura para copiar archivos luego
    T_anterior = T
#############################

