# -*- coding: utf-8 -*-
import os
import re
import sys
import math
import numpy as np
import qsub as qs
import potentials as pots

class Header(object):

    units = "metal"
    boundary = ['p','p','p',]
    dimension = 3
    atom_style = "atomic"

    def __init__(self):
        pass

    def gen_header(self):
        self.commands = "\n\n\n########## HEADER ##########\n"
        self.commands += "units           %s\n"%self.units
        self.commands += "boundary        %s %s %s\n"%(self.boundary[0],self.boundary[1],self.boundary[2])
        self.commands += "dimension        %d\n"%self.dimension
        self.commands += "atom_style      %s\n"%self.atom_style
    
class End_commands(object):

    def __init__(self, np):
        self.np = np

    def gen_commands(self):
        self.commands  = "\n\n\n##########  TAIL  ##########\n"
        self.commands += 'variable natoms equal "count(all)"\n' 
        self.commands +=  'variable teng equal "c_eatoms"\n'
        self.commands +=  'variable length equal "lx/%d"\n'%self.np
        self.commands +=  'variable ecoh equal "v_teng/v_natoms"\n'

        self.commands += "\n\n\n##########  PRINT  ##########\n"
        self.commands +=  "print 'Total energy (eV) = ${teng};'\n"
        self.commands +=  "print 'Number of atoms = ${natoms};'\n"
        self.commands +=  "print 'Lattice constant (Angstoms) = ${length};'\n"
        self.commands +=  "print 'Cohesive energy (eV) = ${ecoh};'\n"

        self.commands += "print 'All done!'\n" 

<<<<<<< HEAD
class Lattice(object):

    def __init__(self, a0, np):
        self.np = np
        self.a0 = a0

    def gen_commands(self):
        self.commands = "\n\n\n########## LATTICE ##########\n"
        self.commands += "lattice "
        self.commands += "        custom %s a1 1.0 0.0 0.0 a2 0.0 1.0 0.0 a3 0.0 0.0 1.0 &\n"%self.a0
        self.commands += "               origin  0. 0. 0. basis 0.0 0.0 0.0 basis 0.5 0.5 0.5\n"

        self.commands += "region "
        self.commands += "        box block 0 %s 0 %s 0 %s\n"%(self.np, self.np, self.np)
        self.commands += "create_box 2 box\n"
        self.commands += "create_atoms    1 region box basis 1 1 basis 2 2\n"

class Thermo(object):

    def __init__(self, logname, s_steps):
        self.logname = logname
        self.s_steps = s_steps

    def gen_commands(self):
        self.commands = "\n\n\n########## THERMO ##########\n"
        self.commands += "log             %s\n"%self.logname
        self.commands += "thermo             %d\n" % self.s_steps
        #self.infile += "thermo_style    custom step atoms pe ke etotal press vol lx ly lz c_eatoms\n"
        self.commands += "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n" 

        # Set options for how thermodynamic information is computed and printed by LAMMPS.
        self.commands += "thermo_modify   line one format float %.16g\n"
=======
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66

class Equil_NPT(object):

    #Paramters of interatomic potential
    pair_style = "meam"
    pair_coeff = "* * meamf U Zr meafile U Zr"
    Temperature = 300
    seed = 12345
    logname = "log"
    infile = ""
   

    def __init__(self, a0, np, infile_name, pair, n_steps = 1000, s_steps = 10, dump = False):
        self.a0 = a0
        self.np = np
        self.infile_name = infile_name
        self.s_steps = s_steps
 
        self.n_steps = n_steps

        self.pair_style = pair["pair_style"]
        self.pair_coeff = pair["pair_coeff"]

        self.header = Header()
        self.end_commands = End_commands(self.np)
<<<<<<< HEAD
        self.lattice = Lattice(self.a0, self.np)
        self.thermo = Thermo (self.logname, self.s_steps)

        self.header.gen_header()
        self.lattice.gen_commands()
        self.thermo.gen_commands()
        

        self.set_potential()

        self.infile += self.header.commands + self.lattice.commands + self.potential
=======
        self.header.gen_header()
        self.gen_lattice()
        self.set_potential()

        self.infile += self.header.commands + self.lattice + self.potential
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66

        #self.infile += "neighbor        2.0 bin\n"
        #self.infile += "neigh_modify    delay 10\n"
       
        self.set_settings(); self.infile += self.setting

        if dump:
            self.infile += "dump            1 all custom 1 dump id type x y z\n"

        self.infile += self.thermo.commands

    def set_temperature (self, t):
        self.Temperature = t

    def set_seed (self, s):
        self.seed = s
      
<<<<<<< HEAD
=======
    def gen_lattice(self):
        self.lattice = "\n\n\n########## LATTICE ##########\n"
        self.lattice += "lattice "
        self.lattice += "        custom %s a1 1.0 0.0 0.0 a2 0.0 1.0 0.0 a3 0.0 0.0 1.0 &\n"%self.a0
        self.lattice += "               origin  0. 0. 0. basis 0.0 0.0 0.0 basis 0.5 0.5 0.5\n"

        self.lattice += "region "
        self.lattice += "        box block 0 %s 0 %s 0 %s\n"%(self.np, self.np, self.np)
        self.lattice += "create_box 2 box\n"
        self.lattice += "create_atoms    1 region box basis 1 1 basis 2 2\n"

>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66
    def set_potential(self):
        self.potential = "\n\n\n########## POTENTIAL ##########\n"
        self.potential += "pair_style      %s\n"%self.pair_style
        self.potential += "pair_coeff      %s\n"%self.pair_coeff

    def set_settings(self):
        self.setting = "\n\n\n########### SETTINGS ############\n"
        self.setting += "compute csym all centro/atom bcc\n" 
        self.setting += "compute eng all pe/atom\n" 
        self.setting += "compute eatoms all reduce sum c_eng\n" 


    def set_fix_npt(self):
        self.fix_npt = "\n\n\n########## FIX ##########\n"
        self.fix_npt += "reset_timestep	0\n"
        self.fix_npt += "timestep 0.001\n"
        self.fix_npt += "velocity all create %f %d mom yes rot no\n"%(self.Temperature, self.seed)
        self.fix_npt += "fix 1 all npt temp %f %f 1 iso 0 0 1 drag .1\n"%(self.Temperature, self.Temperature)


    def set_tail(self):
        self.tail = "\n\nrun %d\n"%self.n_steps
        self.tail += "unfix 1\n"

        self.end_commands.gen_commands()
        self.tail += self.end_commands.commands 

    def gen_infile(self, base_path):# Create input command script for lammps
<<<<<<< HEAD
        self.fname = self.infile_name
=======
        self.fname = self.infile_name + "_equil"
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66

        self.set_fix_npt(); self.infile += self.fix_npt
        self.set_tail();self.infile += self.tail
        fout = open(base_path + "/" + self.fname,'w')
        fout.write(self.infile)

        fout.close()
        return self.fname

class Minimization(object):
    s_steps = 100
    etol = 1.0e-24
    ftol = 1.0e-24
    maxiter =  50000
    maxeval =  100000
    #Paramters of interatomic potential
    pair_style = "meam"
    pair_coeff = "* * meamf U Zr meafile U Zr"
    infile = ""
    relax=0
    box_relax = 0
<<<<<<< HEAD
    logname = "log"
=======
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66

    def __init__(self, a0, np, pair, infile_name, s_steps, dump = False):
        self.a0 = a0
        self.np = np
        self.infile_name = infile_name
        self.s_steps = s_steps

        self.pair_style = pair["pair_style"]
        self.pair_coeff = pair["pair_coeff"]

        self.header = Header()
        self.end_commands = End_commands(self.np)
<<<<<<< HEAD
        self.lattice = Lattice(self.a0, self.np)
        self.thermo = Thermo (self.logname, self.s_steps)

        self.header.gen_header()
        self.lattice.gen_commands()
        self.thermo.gen_commands()

        self.set_potential()

        self.infile += self.header.commands + self.lattice.commands + self.potential
=======

        self.header.gen_header()
        self.gen_lattice()
        self.set_potential()

        self.infile += self.header.commands + self.lattice + self.potential
>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66

        self.infile += "neighbor        2.0 bin\n"
        self.infile += "neigh_modify    delay 10\n"

        self.infile += "# ---------- Define Settings ---------------------\n"
        self.infile += "compute eng all pe/atom\n" 
        self.infile += "compute eatoms all reduce sum c_eng\n" 
<<<<<<< HEAD

        
        if dump:
            self.infile += "dump            1 all custom 1 dump id type x y z\n"

        # genera thermo
        self.infile += self.thermo.commands

    def set_s_steps(self, s):
        self.s_steps = s

    def set_etol(self, s):
        self.etol = s

    def set_ftol(self, s):
        self.ftol = s

    def set_maxiter(self, s):
        self.maxiter = s

    def set_maxeval(self, s):
        self.maxeval = s 

    def set_thermo (self, s_steps):
        self.s_steps = s_steps

    def set_relax (self, relax=0,box_relax=0):
        self.relax=relax
        self.box_relax = box_relax
        
    def set_potential(self):
        self.potential = "\n\n\n########## POTENTIAL ##########\n"
        self.potential += "pair_style      %s\n"%self.pair_style
        self.potential += "pair_coeff      %s\n"%self.pair_coeff

    def set_minimize(self):
        return "minimize %e  %e  %i  %i \n" % (float(self.etol), self.ftol, self.maxiter, self.maxeval)

    def set_tail(self):
        self.end_commands.gen_commands()
        self.tail = self.end_commands.commands 


    def gen_infile(self, path):# Create input command script for lammps
        if self.relax == 1:
            self.fname = path +  "/" + self.infile_name
        else:
            self.fname = path + "/" +  self.infile_name

        self.infile  += "\n\n\n##########  CALCULO ##########\n"
        self.infile  += "reset_timestep 0\n" 
        if self.relax ==0 and self.box_relax ==0:
            self.infile += "run             0\n"
        elif self.relax == 1 and self.box_relax == 0:
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile += self.set_minimize()
        elif self.relax == 1 and self.box_relax ==1:
            self.infile += "fix             1 all box/relax iso 0.0 vmax 0.001\n"
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile +=  self.set_minimize()

        self.set_tail()
        self.infile += self.tail
        fout = open(self.fname,'w')
        fout.write(self.infile)

=======

        self.infile += "\n\n\n########## THERMO ##########\n"
        
        if dump:
            self.infile += "dump            1 all custom 1 dump id type x y z\n"

        self.infile += "log             log\n"
        self.infile += "thermo             %d\n" % self.s_steps
        #self.infile += "thermo_style    custom step atoms pe ke etotal press vol lx ly lz c_eatoms\n"
        self.infile += "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n" 
        self.infile += "thermo_modify   line one format float %.16g\n"


    def set_thermo (self, s_steps):
        self.s_steps = s_steps

    def set_relax (self, relax=0,box_relax=0):
        self.relax=relax
        self.box_relax = box_relax
        
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

    def set_minimize(self):
        return "minimize %e  %e  %i  %i \n" % (float(self.etol), self.ftol, self.maxiter, self.maxeval)

    def set_tail(self):
        self.end_commands.gen_commands()
        self.tail = self.end_commands.commands 


    def gen_infile(self, path):# Create input command script for lammps
        if self.relax == 1:
            self.fname = path +  "/" + self.infile_name + "_relax"
        else:
            self.fname = path + "/" +  self.infile_name + "_static"

        self.infile  += "\n\n\n##########  CALCULO ##########\n"
        self.infile  += "reset_timestep 0\n" 
        if self.relax ==0 and self.box_relax ==0:
            self.infile += "run             0\n"
        elif self.relax == 1 and self.box_relax == 0:
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile += self.set_minimize()
        elif self.relax == 1 and self.box_relax ==1:
            self.infile += "fix             1 all box/relax iso 0.0 vmax 0.001\n"
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile +=  self.set_minimize()

        self.set_tail()
        self.infile += self.tail
        fout = open(self.fname,'w')
        fout.write(self.infile)

>>>>>>> edc0a70ac55bab03af6488b3d8e05a2d93d39f66
        fout.close()
        return self.fname
