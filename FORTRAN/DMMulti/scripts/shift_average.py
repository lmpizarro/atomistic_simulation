#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import sys

archs = ["veloc_fac_1.dat",  "veloc_fac_2.dat",  "veloc_fac_3.dat",  "veloc_ver_1.dat"]


for arch in archs:
    print arch
    with open(arch, 'r') as f:
       first_row = f.readline()
       y = np.loadtxt(f)

    y[0] = y[0] - np.average(y[0])
    y[1] = y[1] - np.average(y[1])
    y[2] = y[2] - np.average(y[2])
    size = np.size(y[:,0])
    f = open(arch, 'w')
    f.write(first_row)
    for i in range (size):
       data = ("%14.8f  %14.8f  %14.8f\n")%(y[i,0], y[i,1], y[i,2],) 
       f.write(data)
    f.close



