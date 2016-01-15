# -*- coding: utf-8 -*-

#
# type: minimization o npt_equil
#
#
dict_config = {"type":"minimization", 
               "qsub_name": "UZrRelax",
               "run_folder": "relax_calc",
               "cant_nodos":8, 
               "infile_name":"infile_relax", 
               "steps":{"save":100, "total":10000},
               "potential":{"type":"meam", "name": "UZr"},
               "lattice":{"param": 3.47, "nper":5}}
