# -*- coding: utf-8 -*-
import os
import sys

root_path = os.path.dirname(os.path.realpath(__file__))
root_libs = os.path.join(root_path, "../libs")
sys.path.append(root_libs)

import gen_runs as gr

#   LEE PARAMETROS DE ENTRADA
# parametros_t.py tiene que existir
import parametros_t as pars


<<<<<<< HEAD
if __name__ == "__main__":
    gr.run_npt(root_path, pars.dict_config)


=======
    run_folder = 'relax_calc' 
    path_run_folder = os.path.join(root_path,run_folder)
    os.makedirs(path_run_folder)

    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs1 = qs.Qsub_sheldon("UZrEquilibration", cant_nodos, "infile_equil")
    qs1.set_run_path (path_run_folder)
    qs1.gen_file()

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  pots.MEAMPOT_UZr("meam", "UZr")
    pot.gen_files(path_run_folder)
    pair = pot.get_lammps_pot()

    # Default lattice parameter
    a0=3.47
    # cantidad de períodos de la red cristalina
    np = 8
    # pasos grabación
    s_steps = 100 
    # cantidad de pasos calculo
    n_steps = 10000
    min = gl.Equil_NPT (a0, np, "infile", pair, n_steps=n_steps, s_steps=s_steps)
    min.set_temperature (100)
    min.set_seed (23456)
    # genera el infile para lammps
    min.gen_infile(path_run_folder)
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66
