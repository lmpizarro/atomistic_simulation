#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import random
import numpy as np

###############################################################################
# ESCRIBE EL ARCHIVO DE ENTRADA DE DINAMICA MOLECULAR 
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
# ESCRIBE EL ARCHIVO 'tabla_temperatura.dat' EN LA CARPETA  temperatura_T
###############################################################################        
def copia_val_medios(T,rho,N,dire,nombre):
    p=[];sp=[];tem=[];U=[]
    # Copio los archivos de la carpeta local y los escribo en un archivo común
    with open('runs_estadistica.dat','r') as valruns:
        valruns.readline() # Primera linea con el encabezado
        for line in valruns:
            col = [float(x) for x in line.split()]
            p.append(col[4])
            sp.append(col[5])
            tem.append(col[6])
            U.append(col[8])
    arch_comun = os.path.join(dire,nombre)
    fa = 1/np.sqrt(N)
    print('Se promediaron {0} de {1} valores para T = {2} | Rho = {3}'.format(len(p),N,T,rho))
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
    valor_str = '     '.join( format(j,'.5E') for j in valor)         
    head_str = '   ' + '         '.join(head)
    if os.path.isfile(arch_comun):
        with open(arch_comun,'a') as comun:
            comun.write( valor_str +'\n')   
    else:
        with open(arch_comun,'a') as comun:
            # Escribe el encabezado si el archivo no existe
            comun.write( head_str + '\n')
            comun.write( valor_str +'\n')   
###############################################################################

###############################################################################
# ESCRIBE EL ARCHIVO CON LOS RESULTADOS DE CADA CORRIDA A UNA DADA TEMPERATURA
# LEE 'val_medios.dat' EN CADA CARPETA RUNXX
# ESCRIBE EL ARCHIVO 'runs_estadistica.dat' EN LA CARPETA T_tmpfolder
# ---- Es usado sólo por corridas.py 
# ---- corridas_paralelo.py usa dos funciones separadas (por la comunicación)
# NO SE MODIFICO PARA DINAMICA MOLECULAR
################################################################################
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
# LEE 'val_medios.dat' EN CADA CARPETA RUNXX 
# LOS ORDENA EN UN STRING Y LOS PASA COMO SALIDA JUNTO AL ENCABEZADO
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
    head = [head1[0]] + [head1[3]] + [head1[1]] + head2 + head3[0:2]
    dato = [dato1[0]] + [dato1[3]] + [dato1[1]] + dato2 + dato3[0:2]
    # Ordena los valores leidos en un string y los pasa como salida
    sep_hea = '          '
    sep_dat = '     '
    head  = 'RUN' + sep_dat + sep_dat.join(head[0:3]) + \
            sep_hea + sep_hea.join(head[3:]) + '\n'
    lin_dat = '%02i'%i + sep_dat + sep_dat.join(dato) + '\n'
    return (lin_dat,head)
###############################################################################

###############################################################################
# ESCRIBE EL ARCHIVO 'runs_estadistica.dat' EN LA CARPETA  temperatura_T
###############################################################################            
def escribe_datos_runs(datos_str,header,path):
    arch_comun = os.path.join(path,'runs_estadistica.dat')
    if os.path.isfile(arch_comun):
        with open(arch_comun,'a') as comun:
            comun.write(datos_str)
    else:
        with open(arch_comun,'a') as comun:
            comun.write(header)
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

###############################################################################
# CONVIERTE DENSIDAD A LADO DE LA CAJA PARA ESCRIBIR AL ARCHIVO DE ENTRADA 
###############################################################################           
def rho_2_lado(n_p,densidad):
    lado = np.power( n_p/densidad , np.float(1)/3 ) 
    lado = str(lado)
    return lado
###############################################################################

