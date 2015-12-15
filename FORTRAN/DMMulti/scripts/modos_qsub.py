#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import errno
import numpy as np
import sys

###############################################################################
# ESCRIBE UN ARCHIVO CON LA SEMILLA PARA ZIGGURAT
###############################################################################           
def escribe_semilla(nombre_arch):
    import random
    random.seed()
    # Obtengo un número aleatorio entero de 4 bytes
    entero_aleatorio = random.randint(-2147483648,2147483647)
    with open(nombre_arch,'w') as fileout:
        fileout.write(str(entero_aleatorio))

cantCorridas = 100000
line1 = "2 1"
line2 = "1.0 1.0 1.0"
line3 = "0.8 0.9 " # 0.6: parte variable
line4 = "5.95 4 2"
line7 = "0.5 0.001" #20000: append variable  
line75 = "10 0.5"
line6 = ".5"

l1qsub = "#PBS -N " # el nombre prueba
qsubc = '''
#PBS -N prueba
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -q verylong
#PBS -o nas-0-0.local:$PBS_O_WORKDIR/$PBS_JOBID.output                                                                       
#PBS -e nas-0-0.local:$PBS_O_WORKDIR/$PBS_JOBID.error                                                                     
cd $PBS_O_WORKDIR
 
cat $PBS_NODEFILE > $PBS_O_WORKDIR/$PBS_JOBID.machines
cat $PBS_NODEFILE | uniq > $PBS_O_WORKDIR/$PBS_JOBID.mpd.hosts

NUM_PROCS=`cat $PBS_O_WORKDIR/$PBS_JOBID.machines|wc -l`
NUM_NODES=`cat $PBS_O_WORKDIR/$PBS_JOBID.mpd.hosts|wc -l`

echo NUM_PROCS = $NUM_PROCS
echo NUM_NODES = $NUM_NODES

export NUM_PROCS 
export NUM_NODES

bin="../dmmulti"
/home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -machinefile $PBS_O_WORKDIR/$PBS_JOBID.machines -np $NUM_PROCS $bin
'''

dirs = ["masa_1.0", "masa_2.0", "masa_3.0", "masa_4.0", "masa_5.0","masa_6.0",\
"masa_7.0", "masa_8.0", "masa_0.6", "masa_1.6", "masa_2.6", "masa_3.6", "masa_4.6","masa_5.6",\
"masa_6.6", "masa_7.6"]

#if size != len(dirs):
#    print "cantidad de nodos inadecuada"
#    sys.exit (0)

# Directorio raíz donde está el ejecutable y este script
curr_dir = os.getcwd()

masas = []
for d in dirs:
    masas.append((d.split("_")[1]))

config = ""

for rank, m in enumerate(masas):
	# lo que corre cada procesador
	print ("rank "), rank, masas[rank]
	config = config +  line1 + "\n"
	config = config + line2 + "\n"
	config = config + line3 + masas[rank] + "\n"
	config = config + line4 + "\n"
	config = config + line7 + " " + str(cantCorridas) + " " + line75 + "\n"
	config = config + line6 + "\n\n"

	dir_corrida = "./" + dirs[rank]

	os.makedirs(dir_corrida)
	f= open(dir_corrida + "/input.par", 'w')
	f.write(config)
	f.close

	# Pasa a carmeta de cada corrida RUN..
	os.chdir(dir_corrida)
	escribe_semilla("seed.dat")
	local_dir = os.getcwd()
        linea1 = l1qsub + "masa_"+ masas[rank]
        qsub = linea1 + qsubc 
	f= open("./run.qsub", 'w')
	f.write(qsub)
	f.close
	os.chdir("../")
        config = ""
