# -*- coding: utf-8 -*-
import os
import re
import sys
import math
import numpy as np

class MEAMPOT_UZr(object):
    meafile ='''rc = 5.5
delta(1,2) = 0.8
re(1,2) = 2.85
delr = 0.1
alpha(1,2) = 4.800000
lattce(1,2) = b1
rho0(1) = 1.2 
rho0(2) = 0.8
ialloy = 1
emb_lin_neg = 1
bkgd_dyn = 1
zbl(1,1) = 1000
zbl(1,2) = 1000
zbl(2,2) =1000
attrac(1,1) = 0.105
repuls(1,1) = 0.105
attrac(1,2) = 0.05
repuls(1,2) = 0.05
attrac(2,2) = -0.03
repuls(2,2) = -0.03
Cmin(1,1,1) = 1.00
Cmax(1,1,1) = 1.9
Cmin(1,1,2) = 0.5
Cmax(1,1,2) = 2.0
Cmin(1,2,1) = 0.8
Cmax(1,2,1) = 2.5
Cmin(1,2,2) = 0.5
Cmax(1,2,2) = 2.8
Cmin(2,2,1) = 0.5
Cmax(2,2,1) = 2.8
Cmin(2,2,2) = 0.7 
Cmax(2,2,2) = 0.99'''

    meamf='''# meam data from vax files fcc,bcc,dia
elt        lat     z       ielement     atwt
alpha      b0      b1      b2           b3    alat    esub    asub
t0         t1              t2           t3            rozero  ibar
'U'        'fcc'    12      92        238.03             
5.1         4.80    6.0     6         6        4.280   5.27   0.98
1           2.50            4         1.0       1       0   
'Zr'       'bcc'    8.      40        91.224
4.10        2.80    2.0     7.00      1        3.535   6.20   0.48
1.0         3.00            2.0      -7       1       0'''

    def __init__(self):
        fout = open('meafile','w')
        fout.write(self.meafile)
        fout.close()

        fout = open('meamf','w')
        fout.write(self.meamf)
        fout.close()


class Qsub(object):
    qsub = '''#PBS -N EquilNPT100
    #PBS -l nodes=8:ppn=2
    #PBS -S /bin/bash
    #PBS -q verylong
    #PBS -o nas-0-0.local:$PBS_O_WORKDIR/$PBS_JOBID.output                                                                   
    #PBS -e nas-0-0.local:$PBS_O_WORKDIR/$PBS_JOBID.error 
    cd $PBS_O_WORKDIR

    cat $PBS_NODEFILE > $PBS_O_WORKDIR/$PBS_JOBID.machines
    cat $PBS_NODEFILE | uniq > $PBS_O_WORKDIR/$PBS_JOBID.mpd.hosts

    NUM_PROCS=`cat $PBS_O_WORKDIR/$PBS_JOBID.machines|wc -l`
    NUM_NODES=`cat $PBS_O_WORKDIR/$PBS_JOBID.mpd.hosts|wc -l`

    echo NUM_PROCS = $NUM_PROCS
    echo NUM_NODES = $NUM_NODES

    export NUM_PROCS 
    export NUM_NODES


    HOME_CALC="/home/pizarro/qsub/python/SimRui/ASE/lammps/UZr/equilibraPy/"
    LOG_FILE="logfile.dat"
    IN_FILE="infile_relax"
    LAMMPS="/home/pizarro/packages/lammps-10Aug15/src/lmp_mpi"

    bin=$LAMMPS" -in "$HOME_CALC$IN_FILE" -l "$HOME_CALC$LOG_FILE
    /home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -machinefile $PBS_O_WORKDIR/$PBS_JOBID.machines -np $NUM_PROCS $bin'''

    qsub_split = qsub.split('\n')

    def __init__(self, name, nodes, infile_name):
        self.name = name
        self.nodes = nodes
        self.infile_name = infile_name

        self.qsub_split[0] ="#PBS -N %s"%self.name  
        self.qsub_split[1] ="#PBS -l nodes=%d:ppn=2"%self.nodes  
        self.qsub_split[21] ="HOME_CALC=%s/"%os.path.dirname(os.path.realpath(__file__))
        self.qsub_split[23] ="IN_FILE=%s"%self.infile_name

        self.file_content = ""
        for e in self.qsub_split:
            self.file_content += e.strip() + "\n"
             

        fout = open('run.qsub','w')
        fout.write(self.file_content)

        fout.close()


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

    def __init__(self, a0, np, infile_name, s_steps, dump = False):
        self.a0 = a0
        self.np = np
        self.infile_name = infile_name
        self.s_steps = s_steps

        self.gen_header()
        self.gen_lattice()
        self.set_potential()

        self.infile += self.header + self.lattice + self.potential

        self.infile += "neighbor        2.0 bin\n"
        self.infile += "neigh_modify    delay 10\n"

        self.infile += "# ---------- Define Settings ---------------------\n"
        self.infile += "compute eng all pe/atom\n" 
        self.infile += "compute eatoms all reduce sum c_eng\n" 

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

        
    def gen_lattice(self):
        self.lattice = "\n\n\n########## LATTICE ##########\n"
        self.lattice += "lattice "
        self.lattice += "        custom %s a1 1.0 0.0 0.0 a2 0.0 1.0 0.0 a3 0.0 0.0 1.0 &\n"%self.a0
        self.lattice += "               origin  0. 0. 0. basis 0.0 0.0 0.0 basis 0.5 0.5 0.5\n"

        self.lattice += "region "
        self.lattice += "        box block 0 %s 0 %s 0 %s\n"%(self.np, self.np, self.np)
        self.lattice += "create_box 2 box\n"
        self.lattice += "create_atoms    1 region box basis 1 1 basis 2 2\n"


    def gen_header(self):
       self.header = "\n\n\n########## HEADER ##########\n"
       self.header += "units           metal\n"
       self.header += "boundary        p p p\n"
       self.header += "dimension        3\n"
       self.header += "atom_style      atomic\n"


    def set_potential(self):
        self.potential = "\n\n\n########## POTENTIAL ##########\n"
        self.potential += "pair_style      %s\n"%self.pair_style
        self.potential += "pair_coeff      %s\n"%self.pair_coeff

    def set_minimize(self):
        return "minimize %e  %e  %i  %i \n" % (float(self.etol), self.ftol, self.maxiter, self.maxeval)


    def set_tail(self):
        self.tail  = "\n\n\n##########  TAIL  ##########\n"
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
        if relax == 1:
            self.fname = self.infile_name + "_relax"
        else:
            self.fname = self.infile_name + "_static"

        self.infile  += "\n\n\n##########  CALCULO ##########\n"
        self.infile  += "reset_timestep 0\n" 
        if relax ==0 and box_relax ==0:
            self.infile += "run             0\n"
        elif relax == 1 and box_relax == 0:
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile += self.set_minimize()
        elif relax == 1 and box_relax ==1:
            self.infile += "fix             1 all box/relax iso 0.0 vmax 0.001\n"
            self.infile += "min_style       cg\n"
            self.infile += "min_modify      line quadratic\n"
            self.infile +=  self.set_minimize()

        self.set_tail()
        self.infile += self.tail
        fout = open(self.fname,'w')
        fout.write(self.infile)

        fout.close()
        return self.fname


if __name__ == "__main__":

    # Default lattice parameter
    a0=3.4862
   
    # cantidad de períodos de la red cristalina
    np = 5
    s_steps = 1000 # pasos grabación
    min = Minimization (a0, np, "infile", s_steps)
    # genera el infile para lammps
    min.gen_infile(1,1)

    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs = Qsub("UZrEquilibration", cant_nodos, "infile_relax")

    # graba los potenciales de acuerdo a como los espera el infile
    pot =  MEAMPOT_UZr()
