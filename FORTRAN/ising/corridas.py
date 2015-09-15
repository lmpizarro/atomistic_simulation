#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import shlex
#from subprocess import call
from subprocess import Popen, PIPE
import errno
import random

# Defino la función para escribir el archivo de datos
def escribe_entrada(nombre,valor):
    with open('parametros.dat','r') as f1:
        for line in f1:
            col = line.split()
            if nombre== 'T':
                new=line.replace(col[2],valor)
            elif nombre == 'N':
                new=line.replace(col[4],valor)
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

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()
# Lista con las temperaturas que se desean calcular
tempe = [0.5, 0.6]
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
    # Se mete en la carpeta
    os.chdir(path_carpeta)       
    # Cambia el archivo de entrada adentro de la carpeta
    escribe_entrada('T',T)
    # Escribe la semmilla en la carpeta
    escribe_semilla()
    # Corre el programa
    #call(curr_dir+'/ising')
    run = Popen(shlex.split(curr_dir+'/ising'),stdout=PIPE)
    # Escribo un archivo log con la info de la corrida
    f=open("log.txt",'w')
    f.write(str(run.communicate()))
    f.close()    
    # Sale de la carpeta
    os.chdir(curr_dir)

