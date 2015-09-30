# -*- coding: utf-8 -*-

#!/usr/bin/python

import os
import shutil
import random
import numpy as np

###############################################################################
# ESCRIBE EL ARCHIVO DE ENTRADA DE ISING
###############################################################################
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
###############################################################################

###############################################################################
# ESCRIBE UN ARCHIVO CON LA SEMILLA PARA ZIGGURAT
###############################################################################           
def escribe_semilla():
    random.seed()
    # Obtengo un número aleatorio entero de 4 bytes
    entero_aleatorio = random.randint(-2147483648,2147483647)
    with open('seed.dat','w') as fileout:
        fileout.write(str(entero_aleatorio))
###############################################################################

###############################################################################
# HACE ESTADISTICA CON LOS RESULTADOS OBTENIDOS A UNA DADA TEMPERATURA
# ESCRIBE EL ARCHIVO 'tablas_temperatura.dat' EN LA CARPETA T_tmpfolder
###############################################################################        
def copia_val_medios(T,N,dire):
    E=[];M=[];c=[];s=[];a=[]
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('runs_estadistica.dat','r') as valruns:
        for line in valruns:
            col = [float(x) for x in line.split()]
            E.append(col[2])
            c.append(col[3])            
            M.append(abs(col[4]))
            s.append(col[5])
            a.append(col[6])
    arch_comun = os.path.join(dire,'tablas_temperatura.dat')
    fa = 1/np.sqrt(N)
    print('Se promediaron {0} de {1} valores para T = {2}'.format(len(E),N,T))
    valor = [float(T), np.mean(E), np.std(E)*fa,
                       np.mean(M), np.std(M)*fa,
                       np.mean(c), np.std(c)*fa,
                       np.mean(s), np.std(s)*fa,
                       np.mean(a), np.std(a)*fa]
    # Abro el archivo donde se volcarán las resultados de cada temperatura   
    with open(arch_comun,'a') as comun:  
        valor_str = '  '.join( format(j,'.5E') for j in valor)         
        comun.write( valor_str +'\n')   
###############################################################################
        
###############################################################################
# ESCRIBE EL ARCHIVO CON LOS RESULTADOS DE CADA CORRIDA A UNA DADA TEMPERATURA
# LEE 'val_medios.dat' EN CADA CARPETA RUNXX
# ESCRIBE EL ARCHIVO 'runs_estadistica.dat' EN LA CARPETA T_tmpfolder
###############################################################################            
def copia_val_medios_runs(i,path):
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('val_medios.dat','r') as vmed:
        valor=vmed.readline()
    # Lee el archivo con la cantidad de aceptaciones
    with open('aceptaciones.dat','r') as archacept:
        acep=archacept.readline()
        colu = [float(x) for x in acep.split()]
        # Guardará el porcentaje de aceptaciones
        por_acep = 100.0*colu[1]/colu[2]
    arch_comun = os.path.join(path,'runs_estadistica.dat')
    with open(arch_comun,'a') as comun:
        comun.write('%02i'%i + ' ' + valor.rstrip() + ' '+ str(por_acep) +'\n') 
###############################################################################

###############################################################################
# BUSCA EL ARCHIVO CON EL ESTADO FINAL A LA TEMPERATURA ANTERIOR Y LO COPIA
# EN LA CARPETA DE LA TEMPERATURA ACTUAL
###############################################################################      
def copia_estado_temp_anterior(path_act,T_ant,T_act):
    # Copia y reescribe el último estado calculado a una temperatura anterior
    if T_ant != []:
        # Reemplazo sólo la parte de la temperatura en toda la ruta
        path_ant = path_act.replace(T_act,T_ant)
        # Ruta completa del archivo a copiar
        # arch_ant = os.path.join(path_ant,'estado.dat')
        # Para todas las corridas utiliza el ultimo estado de RUN00
        arch_ant = os.path.join(os.path.dirname(path_ant),'RUN00','estado.dat')        
        # Copia el archivo
        shutil.copy(arch_ant,path_act)
        # reescribe el archivo
        os.rename('estado.dat','ultimo_estado.dat')
###############################################################################
