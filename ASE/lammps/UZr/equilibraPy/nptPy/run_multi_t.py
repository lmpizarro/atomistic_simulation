#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import subprocess

root_path = os.path.dirname(os.path.realpath(__file__))
root_libs = os.path.join(root_path, "../libs")
sys.path.append(root_libs)

import gen_runs as gr


#   LEE PARAMETROS DE ENTRADA
# parametros_t.py tiene que existir
import parametros_t as pars

info_salida = ""
qsub_present = True
if __name__ == "__main__":
    folders = gr.create_folders_calcs([pars.Temp_min, pars.Temp_max, pars.dTemp])

    # genera el run.qsub para sheldon
    for f in folders:
        # Pasa a carpeta de cada corrida
        os.chdir(f[0])

        if qsub_present:
            #proc = subprocess.Popen(['qsub', './run.qsub'],stdout=subprocess.PIPE)
            proc = subprocess.Popen(['ls', './run.qsub'],stdout=subprocess.PIPE)
            salida = proc.communicate()[0]
            info_salida += salida 
        else:
            print "command to run not configured"


        info_salida += f[0] + "\n"
        os.chdir(root_path)

    fout = open(root_path + "/" + "salida.qsub",'w')
    fout.write(info_salida)
    fout.close()
 
