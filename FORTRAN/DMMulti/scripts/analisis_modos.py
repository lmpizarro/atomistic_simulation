#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

dirs = ["masa_0.6", "masa_1.6", "masa_2.6", "masa_3.6", "masa_4.6","masa_5.6",\
"masa_6.6", "masa_7.6", "masa_8.6"]

names = ["modo_fac_1.dat", "modo_fac_2.dat", "modo_fac_3.dat", "modo_ver_1.dat"]


for n in names:
    print n
    for d in dirs:
        arch = "./" + d + "/" + n
        y = np.loadtxt(arch,skiprows=1)
        salida = np.array([np.average(y[0]), np.average( y[1]), np.average(y[2])])
        print("%s %14.8f %14.8f %14.8f %14.8f" )%(d.split("_")[1], salida[0], salida[1],
                salida[2], np.average(salida))
    print ("")    


