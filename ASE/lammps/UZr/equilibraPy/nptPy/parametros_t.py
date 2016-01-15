# -*- coding: utf-8 -*-

###############################################################################       
#   PARAMETROS DE ENTRADA PARA CORRER A DISTINTAS TEMPERATURAS
###############################################################################
#------ Barrido de temperaturas
# Temperatura mínima
Temp_min = 1100 
# Temperatura máxima
Temp_max = 1400
# Paso de temperatura
dTemp = 100 
# Default lattice parameter

#
# type: minimization o npt_equil
#
#
dict_config = {"type":"npt_equil", 
               "qsub_name": "UZrEquil",
               "run_folder": "calc_npt",
               "cant_nodos":8, 
               "infile_name":"infile_equil", 
               "steps":{"save":100, "total":10000},
               "Temperature": 1200,
               "seed": 45678,
               "potential":{"type":"meam", "name": "UZr"},
               "lattice":{"param": 3.47, "nper":5}}
