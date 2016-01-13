#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import re
import sys
import math
import numpy as np
import errno
import subprocess

root_path = os.path.dirname(os.path.realpath(__file__))
root_libs = os.path.join(root_path, "../libs")
sys.path.append(root_libs)

import qsub as qs
import potentials as pots
import gen_lammps as gl


#   LEE PARAMETROS DE ENTRADA
if os.path.isfile("parametros_t.py"):
    import parametros_t as parametros

    Temp_min = np.float(parametros.Temp_min)
    Temp_max = np.float(parametros.Temp_max)
    dTemp    = np.float(parametros.dTemp)
    a0 =  np.float(parametros.a0)
    # cantidad de períodos de la red cristalina
    nper = int(parametros.nper) 
    # pasos grabación
    s_steps =  int(parametros.s_steps)
    # cantidad de pasos calculo
    n_steps = int(parametros.n_steps)


else:
    Temp_min = 1000 
    Temp_max = 1400
    dTemp = 50 
    # Default lattice parameter
    a0=3.47
    # cantidad de períodos de la red cristalina
    nper = 4
    # pasos grabación
    s_steps = 100 
    # cantidad de pasos calculo
    n_steps = 10000

folders = []
def create_folders_calcs(tempes):
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
        # Se crea la carpeta sólo si ésta no existe. De lo contrario se saltea la 
        # temperatura correspondiente
        try: 
            os.makedirs(path_carpeta)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            else:
                print('Ya existe el directorio. Se omite corrida con T=' + Tnombre)
                T_anterior = Tnombre
                # Si se llega al final de la lista, se finaliza el programa
                if T==tempe[-1]:
                    print('No hay nada por hacer, se sale del programa')
                    sys.exit()
                continue

info_salida = ""
if __name__ == "__main__":
    create_folders_calcs([Temp_min, Temp_max, dTemp])


    pot =  pots.MEAMPOT_UZr("meam", "UZr")
    pair = pot.get_lammps_pot()



    # genera el run.qsub para sheldon
    cant_nodos = 8
    for f in folders:
        qs1 = qs.Qsub_sheldon("UZrEquil" + str(f[1]), cant_nodos, "infile_equil")
        qs1.set_run_path (f[0])
        qs1.gen_file()
        # graba los potenciales de acuerdo a como los espera el infile
        pot.gen_files(f[0])

        min = gl.Equil_NPT (a0, nper, "infile", pair, n_steps=n_steps, s_steps=s_steps)
        min.set_seed (23456)

        print f[0], f[1]
        min.set_temperature (f[1])
        min.gen_infile(f[0], 1, 1)

        # Pasa a carmeta de cada corrida RUN..
        os.chdir(f[0])

        proc = subprocess.Popen(['qsub', './run.qsub'],stdout=subprocess.PIPE)
        salida = proc.communicate()[0]
        info_salida += salida 
        os.chdir(root_path)

    fout = open(root_path + "/" + "salida.qsub",'w')
    fout.write(info_salida)
    fout.close()
 
