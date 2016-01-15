# -*- coding: utf-8 -*-

import os

import qsub as qs
import potentials as pots
import gen_lammps as gl
import numpy as np

def gen_mpi (path, nper, infile):
    command = "mpirun -np %d lmp_mpi -in %s" % (nper, infile)
    fout = open(path + '/run.sh','w')
    fout.write(command)
    fout.close()

def create_folders_calcs(tempes):
    folders = []
    Temp_min = tempes[0]
    Temp_max = tempes[1]
    dTemp = tempes[2]
    # Lista de temperaturas
    tempe = np.arange(Temp_min, Temp_max + dTemp, dTemp)
    # ordena de mayor a menor el array numpy
    tempe = np.fliplr([tempe])[0]

    # Directorio raíz donde está el ejecutable y este script
    curr_dir = os.getcwd()

    for T in tempe:
        Tnombre = str(T)
        # Nombre de la carpeta uqe se va a crear
        carpeta = 'temperatura_' + Tnombre 
        # Camino completo de la carpeta que se va a crear
        path_carpeta = os.path.join(curr_dir,carpeta)
        folders.append([path_carpeta,T])
    return folders


def run_npt(root_path, dict_config):

    run_folder = dict_config["run_folder"] 
    path_run_folder = os.path.join(root_path,run_folder)
    os.makedirs(path_run_folder)

    # genera el run.qsub para sheldon
    cant_nodos = dict_config["cant_nodos"]
    qs1 = qs.Qsub_sheldon( dict_config["qsub_name"], 
                          cant_nodos, 
                          dict_config["infile_name"])
    qs1.set_run_path (path_run_folder)
    qs1.gen_file()

    gen_mpi(run_folder, dict_config["cant_nodos"], dict_config["infile_name"])

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  pots.MEAMPOT_UZr(dict_config["potential"]["type"], 
                            dict_config["potential"]["name"])
    pot.gen_files(path_run_folder)
    pair = pot.get_lammps_pot()

    # Default lattice parameter
    a0=dict_config["lattice"]["param"]
    # cantidad de períodos de la red cristalina
    nper = dict_config["lattice"]["nper"]
    # pasos grabación
    s_steps = dict_config["steps"]["save"] 
    # cantidad de pasos calculo
    n_steps =  dict_config["steps"]["total"]
    npt = gl.Equil_NPT (a0, nper, dict_config["infile_name"], pair, n_steps=n_steps, s_steps=s_steps)
    npt.set_temperature (dict_config["Temperature"])
    npt.set_seed (dict_config["seed"])
    # genera el infile para lammps
    npt.gen_infile(path_run_folder)


def run_relax(root_path, dict_config):
    run_folder = dict_config["run_folder"] 
    path_run_folder = os.path.join(root_path,run_folder)
    os.makedirs(path_run_folder)

    # genera el run.qsub para sheldon
    qs1 = qs.Qsub_sheldon(dict_config["qsub_name"],
                          dict_config["cant_nodos"], 
                          dict_config["infile_name"])
    qs1.set_run_path (path_run_folder)
    qs1.gen_file()

    gen_mpi(run_folder, dict_config["cant_nodos"], dict_config["infile_name"])

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  pots.MEAMPOT_UZr(dict_config["potential"]["type"], 
                            dict_config["potential"]["name"])
    pot.gen_files(path_run_folder)
    pair = pot.get_lammps_pot()

    # Default lattice parameter
    a0= dict_config["lattice"]["param"]
    # cantidad de períodos de la red cristalina
    nper =  dict_config["lattice"]["nper"]

    s_steps =  dict_config["steps"]["save"] # pasos grabación
    min = gl.Minimization (a0, nper, pair, dict_config["infile_name"], s_steps)
    min.set_relax(1,1)
    # genera el infile para lammps
    min.gen_infile(path_run_folder)


