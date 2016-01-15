# -*- coding: utf-8 -*-
import os
import sys


root_path = os.path.dirname(os.path.realpath(__file__))
root_libs = os.path.join(root_path, "../libs")
sys.path.append(root_libs)

import qsub as qs
import potentials as pots
import gen_lammps as gl
import gen_runs as gr

import parametros_t as pars


if __name__ == "__main__":
    gr.run_relax(root_path, pars.dict_config)


