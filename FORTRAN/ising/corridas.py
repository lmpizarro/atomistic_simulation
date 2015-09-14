#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
from subprocess import call
import errno

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()

tempe = [0.5, 0.6]
tempe = [str(i) for i in tempe]

# Loop para crear todos los directorios y correr el ejecutable en ellos
for T in tempe:
    with open('parametros.dat','r') as f1:
        for line in f1:
            col = line.split()
            new=line.replace(col[2],T)
   
    with open('parametros.dat','w') as fileout:
        fileout.write(new)

    carpeta = T + '_tmpfolder'

    path_carpeta = os.path.join(curr_dir,carpeta)

    try: 
        os.makedirs(path_carpeta)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print('Ya existe el directorio. Se omite corrida con T=' + T)
            continue
        
    shutil.copy('parametros.dat',path_carpeta) 
    os.chdir(path_carpeta)       
    call(curr_dir+'/ising')
    
    os.chdir(curr_dir)
    # Creo directorio
    # copio ejecutable y entradas
    # Hago algo con la salida
    