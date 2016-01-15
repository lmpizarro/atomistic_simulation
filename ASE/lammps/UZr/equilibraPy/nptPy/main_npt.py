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


if __name__ == "__main__":
    gr.run_npt(root_path, pars.dict_config)


