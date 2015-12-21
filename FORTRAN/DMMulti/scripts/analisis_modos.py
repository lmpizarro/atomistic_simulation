#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

dirs = ["masa_0.6", "masa_1.0", "masa_1.6", "masa_2.0", "masa_2.6", "masa_3.0",\
        "masa_3.6", "masa_4.0", "masa_4.6", "masa_5.0", "masa_5.6", "masa_6.0",\
        "masa_6.6", "masa_7.0", "masa_7.6", "masa_8.0"]

r_name = "veloc"
names = ["_fac_1.dat", "_fac_2.dat", "_fac_3.dat", "_ver_1.dat"]


for n in names:
    print r_name + n
    for d in dirs:
        arch = "./" + d + "/" + r_name + n
        y = np.loadtxt(arch,skiprows=1)
        salida = np.array([np.average(y[0]), np.average( y[1]), np.average(y[2])])
        print("%s %14.8f %14.8f %14.8f %14.8f" )%(d.split("_")[1], salida[0], salida[1],
                salida[2], np.average(salida))
    print ("")    


