#!/usr/bin/python
# -*- coding: utf-8 -*-

# Crea una carpeta para cada temperatura. En cada una de esas carpetas, a su
# vez, crea Nruns carpetas para hacer estadística y obtener los valores con
# sus respectivos errores.

import os
import shutil
import subprocess
import errno
import random
import numpy as np

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
    valor = [float(T), np.mean(E), np.std(E)*fa, np.mean(M), np.std(M)*fa,
                np.mean(c), np.std(c)*fa, np.mean(s), np.std(s)*fa,
                np.mean(a), np.std(a)]
    # Abro el archivo donde se volcarán las resultados de cada temperatura   
    with open(arch_comun,'a') as comun:  
        valor_str = '  '.join( ('{:.5E}').format(j) for j in valor) 
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

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()
# Lista con las temperaturas que se desean calcular
T_min = 0.5
T_max = 1.0
dT    = 0.1
tempe = np.arange(T_min,T_max+dT,dT)
#tempe = tempe.tolist() + [2.15, 2.25, 2.35, 2.45, 2.55] # Mäs detalle en la Tc
tempe.sort()
# lo convierto a string
tempe = [str(i) for i in tempe]
# Número de corridas para cada temperatura
Nrun  = 3

# Loop para crear todos los directorios y correr el ejecutable en ellos
for T in tempe:
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
    # Corre el programa para ver la convergencia
    escribe_entrada('N','4000') 

   # Loop para correr N veces con los mismos parámetros para calcular el error
    for i in range(0,Nrun):
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
        salida = subprocess.check_output(curr_dir+'/ising')
    
        # Guardo la salida para ver ue hizo
        f=open("log1.txt",'w')
        f.write(salida)
        f.close() 
        # Guardo las salidas por si hacen falta
        os.rename('energia.dat','energia_terma.dat')
        os.rename('magneti.dat','magneti_terma.dat')    
    
        # Corre por segunda vez tomando el estado anterior. Aumento el N
        escribe_entrada('N','10000')
        print('Corriendo {} a la temperatura {}'.format(carpeta_runs,T))
        salida = subprocess.check_output(curr_dir+'/ising')
    
        # Guardo la salida para ver ue hizo
        f=open("log2.txt",'w')
        f.write(salida)
        f.close()    
        # Copia los valores medios
        copia_val_medios_runs(i)
        # Sale de la carpeta de corridas        
        os.chdir(path_carpeta)
    # Hace estadística de todas las corridas y lo copia en un archivo
    copia_val_medios(T,Nrun)
    # Sale de la carpeta
    os.chdir(curr_dir)

#############################

