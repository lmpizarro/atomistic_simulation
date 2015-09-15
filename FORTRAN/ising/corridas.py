#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import shlex
import subprocess
import errno
import random

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
            

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()
# Lista con las temperaturas que se desean calcular
tempe = [x*0.1 for x in range(5,50)]
tempe = tempe + [2.15, 2.25, 2.35, 2.45, 2.55]

tempe.sort()

tempe = [str(i) for i in tempe]

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
    # Escribe la semmilla en la carpeta
    escribe_semilla()

    # Corre el programa para ver la convergencia
    escribe_entrada('N','100000') 
    salida = subprocess.check_output(curr_dir+'/ising')
    
    # Guardo la salida para ver ue hizo
    f=open("log1.txt",'w')
    f.write(salida)
    f.close() 
    
    # Corre por segunda veztomando el estado anterior. Aumento el N
    escribe_entrada('N','3000000')
    salida = subprocess.check_output(curr_dir+'/ising')
    
    # Guardo la salida para ver ue hizo
    f=open("log2.txt",'w')
    f.write(salida)
    f.close()    

    # Copia los valores medios
    copia_val_medios()
    # Sale de la carpeta
    os.chdir(curr_dir)

#############################

