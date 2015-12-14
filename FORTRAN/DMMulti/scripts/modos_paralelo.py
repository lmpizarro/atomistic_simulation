#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import errno
import numpy as np
import sys

# Paralelización
from mpi4py import MPI
#from mpi4py.MPI import ANY_SOURCE

# Información de dónde está            
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

line1 = "2 1"
line2 = "1.0 1.0 1.0"
line3 = "0.8 0.9 " # 0.6: parte variable
line4 = "5.95 4 2"
line7 = "0.5 0.001" #20000: append variable  
line75 = "10 0.5"
line6 = ".5"

dirs = ["masa_0.6", "masa_1.6", "masa_2.6", "masa_3.6", "masa_4.6","masa_5.6",\
"masa_6.6", "masa_7.6"]

if size != len(dirs):
    print "cantidad de nodos inadecuada"
    sys.exit (0)

masas = []
for d in dirs:
    masas.append((d.split("_")[1]))



if rank==0:
    print ("rank 0")
# Todos esperan a que root haya terminado
comm.Barrier()

# lo que corre cada procesador
print ("rank "), rank, dir[rank]
print line1
print line2
print line3 + str(rank)
print line4
print line7 + " 20000 " + line75
print line6

