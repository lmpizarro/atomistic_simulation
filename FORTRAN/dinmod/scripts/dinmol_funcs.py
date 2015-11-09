#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import random
import numpy as np

###############################################################################
# ESCRIBE EL ARCHIVO DE ENTRADA DE ISING
###############################################################################
def escribe_entrada(nombre,valor):
    with open('parametros.dat','r') as f1:
        fila  = f1.readline()
        resto = f1.read()
        col   = fila.split()
        if nombre== 'Temp':
            col[0] = valor
        elif nombre == 'N_part':
            col[1] = valor
        elif nombre == 'L':
            col[2] = valor
        elif nombre == 'dt':
            col[3] = valor
        elif nombre == 'N_pasos':
            col[4] = valor
        elif nombre == 'N_grab':
            col[8] = valor
        else:
            print('No se reconoce el valor a escribir')
            exit
    new = ' '.join(col)
    with open('parametros.dat','w') as fileout:
        fileout.write(new+'\n')
        fileout.write(resto)
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
def copia_val_medios(T,rho,N,dire):
    p=[];sp=[];tem=[];U=[]
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('runs_estadistica.dat','r') as valruns:
        for line in valruns:
            col = [float(x) for x in line.split()]
            p.append(col[4])
            sp.append(col[5])            
            tem.append(col[6])
            U.append(col[8])
    arch_comun = os.path.join(dire,'tablas_temperatura.dat')
    fa = 1/np.sqrt(N)
    print('Se promediaron {0} de {1} valores para T = {2}'.format(len(p),N,T))
    valor = [float(T), float(rho),  
             np.mean(p), np.std(p)*fa,
             np.mean(sp), np.std(sp)*fa,
             np.mean(tem), np.std(tem)*fa,
             np.mean(U), np.std(U)*fa ]
    # Títulos de las columnas escritas
    head = ['Temp', 'Densidad', 
            'mean(p)', 'std(p)',
            'mean(sp)', 'std(sp)',
            'mean(temp)', 'std(temp)',
            'mean(U)', 'std(U)'] 
    # Abro el archivo donde se volcarán las resultados de cada temperatura   
    with open(arch_comun,'a') as comun:  
        #head_str = '       '.join(head)
        #comun.write( head_str + '\n')
        valor_str = '  '.join( format(j,'.5E') for j in valor)         
        comun.write( valor_str +'\n')   
###############################################################################

###############################################################################
# ESCRIBE EL ARCHIVO CON LOS RESULTADOS DE CADA CORRIDA A UNA DADA TEMPERATURA
# LEE 'val_medios.dat' EN CADA CARPETA RUNXX
# ESCRIBE EL ARCHIVO 'runs_estadistica.dat' EN LA CARPETA T_tmpfolder
# ---- Es usado sólo por corridas.py 
# ---- corridas_paralelo.py usa dos funciones separadas (por la comunicación)
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
# LEE 'val_medios.dat' Y 'aceptaciones.dat' EN CADA CARPETA RUNXX 
# LOS ORDENA EN UN STRING Y LOS PASA COMO SALIDA
###############################################################################            

def lee_datos_runs(i):
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('val_medios.dat','r') as vmed:
        head1 = [x for x in vmed.readline().split()]
        dato1 = [x for x in vmed.readline().split()]
        vmed.readline()
        head2 = [x for x in vmed.readline().split()]
        dato2 = [x for x in vmed.readline().split()]
        vmed.readline()
        head3 = [x for x in vmed.readline().split()]
        dato3 = [x for x in vmed.readline().split()]
    # Se arma el string con encabezado y datos de interés para
    # ser escritos en el archivo runs_estadistica.dat
    #head = [head1[0]] + [head1[5]] + head2 + head3[0:2]
    dato = [dato1[0]] + [dato1[3]] + [dato1[1]] + dato2 + dato3[0:2]
    # Ordena los valores leidos en un string y los pasa como salida
    sep = '     '
    #lin_dat ='RUN' + sep + sep.join(head) + '\n'
    lin_dat = '%02i'%i + sep + sep.join(dato) + '\n'
    return lin_dat
###############################################################################

###############################################################################
# ESCRIBE EL ARCHIVO 'runs_estadistica.dat' EN LA CARPETA T_tmpfolder
###############################################################################            

def escribe_datos_runs(datos_str,path):
    arch_comun = os.path.join(path,'runs_estadistica.dat')
    with open(arch_comun,'a') as comun:
        comun.write(datos_str)
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
        # Para todas las corridas utiliza el ultimo estado de RUN00
        arch_ant = os.path.join(path_ant,'RUN00','estados.dat')        
        # Copia el archivo
        shutil.copy(arch_ant,path_act)
###############################################################################