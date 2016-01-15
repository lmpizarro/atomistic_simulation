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


folders = []
if __name__ == "__main__":
    folders = gr.create_folders_calcs([pars.Temp_min, pars.Temp_max, pars.dTemp])

    for f in folders:
        pars.dict_config["run_folder"] = f[0]
        pars.dict_config["Temperature"] = f[1]
        gr.run_npt(root_path, pars.dict_config)
