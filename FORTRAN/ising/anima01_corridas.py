#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess
import shutil

N_tot = 1000
N_tol = 100
direc = 'anima'
j=0
aux = 0
mag = []
mag_anterior = []
os.makedirs(direc)
argu = ['-k'] + ['1'] + ['-t'] + ['0.7']
while (j<N_tot):
    proc = subprocess.Popen([os.getcwd() + '/ising'] + argu,stdout=subprocess.PIPE)
    proc.wait()
    with open('val_medios.dat') as arch:
        line = arch.readline()
        cols = line.split()
    if (cols[3] != mag_anterior) and (mag_anterior!=[]):
        print(cols[3])
        mag.append(cols[3])        
        j += 1
        shutil.copy('estado_mat.dat','./'+direc+'/estado%05i'%j)
        os.rename('estado.dat','ultimo_estado.dat')
        mag_anterior = cols[3]
        aux = 0
    else:
        mag_anterior = cols[3]
        aux += 1
        if aux > N_tol:
            print('No cambia luego de {0} iteraciones'.format(N_tol))
            N_tot = j
            break

# Guardo los datps para poder graficar luego
with open('anim_datos.dat','w+') as arch:
    arch.write(str(N_tot)+'\n')
    for num in mag:
        arch.write(str(num) +'\n')
