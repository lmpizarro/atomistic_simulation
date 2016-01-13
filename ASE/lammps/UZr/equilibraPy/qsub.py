# -*- coding: utf-8 -*-

import os

class Qsub(object):
    run_path = os.path.dirname(os.path.realpath(__file__))
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
        self.qsub_split[21] ="HOME_CALC=%s/"%self.run_path
        self.qsub_split[23] ="IN_FILE=%s"%self.infile_name

    def set_run_path(self, p):
        self.run_path = p 
        self.qsub_split[21] ="HOME_CALC=%s/"%self.run_path

    def gen_file(self):
        self.file_content = ""
        for e in self.qsub_split:
            self.file_content += e.strip() + "\n"
             

        fout = open('run.qsub','w')
        fout.write(self.file_content)

        fout.close()


if __name__ == "__main__":
    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs = Qsub("UZrEquilibration", cant_nodos, "infile_equil")

    run_path = os.path.dirname(os.path.realpath(__file__))
    qs.set_run_path (run_path)
    qs.gen_file()


