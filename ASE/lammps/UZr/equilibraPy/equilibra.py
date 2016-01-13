# -*- coding: utf-8 -*-
import os
import re
import sys
import math
import numpy as np
import qsub as qs
import potentials as pots

class Equilibra(object):

    #Paramters of interatomic potential
    pair_style = "meam"
    pair_coeff = "* * meamf U Zr meafile U Zr"
    Temperature = 300
    seed = 12345
    infile = ""
   

    def __init__(self, a0, np, infile_name, n_steps = 1000, s_steps = 10, dump = False):
        self.a0 = a0
        self.np = np
        self.infile_name = infile_name
        self.s_steps = s_steps
 
        self.n_steps = n_steps

        self.gen_header()
        self.gen_lattice()
        self.set_potential()

        self.infile += self.header + self.lattice + self.potential

        #self.infile += "neighbor        2.0 bin\n"
        #self.infile += "neigh_modify    delay 10\n"
       
        self.set_settings(); self.infile += self.setting

        if dump:
            self.infile += "dump            1 all custom 1 dump id type x y z\n"

        self.set_thermo(); self.infile += self.thermo

    def set_temperature (self, t):
        self.Temperature = t

    def set_seed (self, s):
        self.seed = s

    def gen_header(self):
       self.header = "\n\n\n########## HEADER ##########\n"
       self.header += "units           metal\n"
       self.header += "boundary        p p p\n"
       self.header += "dimension        3\n"
       self.header += "atom_style      atomic\n"
       
    def gen_lattice(self):
        self.lattice = "\n\n\n########## LATTICE ##########\n"
        self.lattice += "lattice "
        self.lattice += "        custom %s a1 1.0 0.0 0.0 a2 0.0 1.0 0.0 a3 0.0 0.0 1.0 &\n"%self.a0
        self.lattice += "               origin  0. 0. 0. basis 0.0 0.0 0.0 basis 0.5 0.5 0.5\n"

        self.lattice += "region "
        self.lattice += "        box block 0 %s 0 %s 0 %s\n"%(self.np, self.np, self.np)
        self.lattice += "create_box 2 box\n"
        self.lattice += "create_atoms    1 region box basis 1 1 basis 2 2\n"

    def set_potential(self):
        self.potential = "\n\n\n########## POTENTIAL ##########\n"
        self.potential += "pair_style      %s\n"%self.pair_style
        self.potential += "pair_coeff      %s\n"%self.pair_coeff

    def set_settings(self):
        self.setting = "\n\n\n########### SETTINGS ############\n"
        self.setting += "compute csym all centro/atom bcc\n" 
        self.setting += "compute eng all pe/atom\n" 
        self.setting += "compute eatoms all reduce sum c_eng\n" 

    def set_thermo (self):
        self.thermo = "\n\n\n########## THERMO ##########\n"
        self.thermo += "log             log\n"
        self.thermo += "thermo             %d\n" % self.s_steps
        #self.infile += "thermo_style    custom step atoms pe ke etotal press vol lx ly lz c_eatoms\n"
        self.thermo += "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n" 

        # Set options for how thermodynamic information is computed and printed by LAMMPS.
        self.thermo += "thermo_modify   line one format float %.16g\n"

    def set_fix_npt(self):
        self.fix_npt = "\n\n\n########## FIX ##########\n"
        self.fix_npt += "reset_timestep	0\n"
        self.fix_npt += "timestep 0.001\n"
        self.fix_npt += "velocity all create %f %d mom yes rot no\n"%(self.Temperature, self.seed)
        self.fix_npt += "fix 1 all npt temp %f %f 1 iso 0 0 1 drag .1\n"%(self.Temperature, self.Temperature)



    # viene de minimiza ver si queda para generalizar todo
    def set_minimize(self):
        return "minimize %e  %e  %i  %i \n" % (float(self.etol), self.ftol, self.maxiter, self.maxeval)


    def set_tail(self):
        self.tail  = "\n\n\n##########  TAIL  ##########\n"
        self.tail += "run %d\n"%self.n_steps
        self.tail += "unfix 1\n"
        self.tail += 'variable natoms equal "count(all)"\n' 
        self.tail +=  'variable teng equal "c_eatoms"\n'
        self.tail +=  'variable length equal "lx/%d"\n'%self.np
        self.tail +=  'variable ecoh equal "v_teng/v_natoms"\n'

        self.tail += "\n\n\n##########  PRINT  ##########\n"
        self.tail +=  "print 'Total energy (eV) = ${teng};'\n"
        self.tail +=  "print 'Number of atoms = ${natoms};'\n"
        self.tail +=  "print 'Lattice constant (Angstoms) = ${length};'\n"
        self.tail +=  "print 'Cohesive energy (eV) = ${ecoh};'\n"

        self.tail += "print 'All done!'\n" 

    def gen_infile(self, relax=0,box_relax=0):# Create input command script for lammps
        self.fname = self.infile_name + "_equil"

        self.set_fix_npt(); self.infile += self.fix_npt
        self.set_tail();self.infile += self.tail
        fout = open(self.fname,'w')
        fout.write(self.infile)

        fout.close()
        return self.fname


if __name__ == "__main__":

    # Default lattice parameter
    a0=3.4862
   
    # cantidad de períodos de la red cristalina
    np = 8
    s_steps = 100 # pasos grabación
    n_steps = 10000
    min = Equilibra (a0, np, "infile", n_steps=n_steps, s_steps=s_steps)
    min.set_temperature (1200)
    min.set_seed (23456)
    # genera el infile para lammps
    min.gen_infile(1,1)

    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs1 = qs.Qsub("UZrEquilibration", cant_nodos, "infile_equil")
    run_path = os.path.dirname(os.path.realpath(__file__))
    qs1.set_run_path (run_path)
    qs1.gen_file()

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  pots.MEAMPOT_UZr("meam", "UZr")
    pot.gen_files()
