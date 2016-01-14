# -*- coding: utf-8 -*-
import os
import re
import sys
import math
import numpy as np


root_path = os.path.dirname(os.path.realpath(__file__))
root_libs = os.path.join(root_path, "../libs")
sys.path.append(root_libs)

import qsub as qs
import potentials as pots
import gen_lammps as gl


if __name__ == "__main__":

    run_folder = 'relax_calc' 
    path_run_folder = os.path.join(root_path,run_folder)
    os.makedirs(path_run_folder)

    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs1 = qs.Qsub_sheldon("UZrRelax", cant_nodos, "infile_relax")
    qs1.set_run_path (path_run_folder)
    qs1.gen_file()

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  pots.MEAMPOT_UZr("meam", "UZr")
    pot.gen_files(path_run_folder)
    pair = pot.get_lammps_pot()

    # Default lattice parameter
    a0=3.4862
   
    # cantidad de períodos de la red cristalina
    np = 5
    s_steps = 1000 # pasos grabación
    min = gl.Minimization (a0, np, pair, "infile", s_steps)
    min.set_relax(1,1)
    # genera el infile para lammps
    min.gen_infile(path_run_folder)


