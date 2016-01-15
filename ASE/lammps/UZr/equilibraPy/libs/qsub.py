# -*- coding: utf-8 -*-

import os

class Qsub_sheldon(object):
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
             

        fout = open(self.run_path + '/run.qsub','w')
        fout.write(self.file_content)

        fout.close()

class Qsub_isaac(object):
    run_path = os.path.dirname(os.path.realpath(__file__))
    qsub =''' ### PONER el run.qsub de referencia### '''

    def __init__(self, name, nodes, infile_name):
        pass

    def set_run_path(self, p):
        pass
 
    def gen_file(self):
        pass

# Para conocer el estado de la cola
# de procesos en sheldon
''' PEXPECT LICENSE

    This license is approved by the OSI and FSF as GPL-compatible.
        http://opensource.org/licenses/isc-license.txt

    Copyright (c) 2012, Noah Spurrier <noah@noah.org>
    PERMISSION TO USE, COPY, MODIFY, AND/OR DISTRIBUTE THIS SOFTWARE FOR ANY
    PURPOSE WITH OR WITHOUT FEE IS HEREBY GRANTED, PROVIDED THAT THE ABOVE
    COPYRIGHT NOTICE AND THIS PERMISSION NOTICE APPEAR IN ALL COPIES.
    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
    WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
    ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
    WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
    ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
    OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

'''

import os, sys, re, getopt, getpass
import pexpect



#
# Some constants.
#
COMMAND_PROMPT = '[#$] ' ### This is way too simple for industrial use -- we will change is ASAP.
TERMINAL_PROMPT = '(?i)terminal type\?'
TERMINAL_TYPE = 'vt100'
# This is the prompt we get if SSH does not have the remote host's public key stored in the cache.
SSH_NEWKEY = '(?i)are you sure you want to continue connecting'


def qstat_u(host_access, commands):
    global COMMAND_PROMPT, TERMINAL_PROMPT, TERMINAL_TYPE, SSH_NEWKEY, HOST,\
    USER, PASSWORD
    ######################################################################
    ## Parse the options, arguments, get ready, etc.
    ######################################################################

    #
    # Login via SSH
    #
    user = host_access["user"]
    host = host_access["host"]
    password = host_access["password"]

    child = pexpect.spawn('ssh -l %s %s'%(user, host))
    i = child.expect([pexpect.TIMEOUT, SSH_NEWKEY, COMMAND_PROMPT, '(?i)password'])
    if i == 0: # Timeout
        error = 'ERROR! could not login with SSH. Here is what SSH said:\n'
        #print(child.before, child.after)
        #print(str(child))
        commands["error": error]
        return commands
    if i == 1: # In this case SSH does not have the public key cached.
        child.sendline ('yes')
        child.expect ('(?i)password')
    if i == 2:
        # This may happen if a public key was setup to automatically login.
        # But beware, the COMMAND_PROMPT at this point is very trivial and
        # could be fooled by some output in the MOTD or login message.
        pass
    if i == 3:
        child.sendline(password)
        # Now we are either at the command prompt or
        # the login process is asking for our terminal type.
        i = child.expect ([COMMAND_PROMPT, TERMINAL_PROMPT])
        if i == 1:
            child.sendline (TERMINAL_TYPE)
            child.expect (COMMAND_PROMPT)
    #
    # Set command prompt to something more unique.
    #
    COMMAND_PROMPT = "\[PEXPECT\]\$ "
    child.sendline ("PS1='[PEXPECT]\$ '") # In case of sh-style
    i = child.expect ([pexpect.TIMEOUT, COMMAND_PROMPT], timeout=10)
    if i == 0:
        #print("# Couldn't set sh-style prompt -- trying csh-style.")
        child.sendline ("set prompt='[PEXPECT]\$ '")
        i = child.expect ([pexpect.TIMEOUT, COMMAND_PROMPT], timeout=10)
        if i == 0:
            error = "Failed to set command prompt using sh or csh style.\n"
            error += "Response was:\n"
            error += child.before + "\n"
            commands["error": error]
            return commands

    for command in commands:
        if command["desc"] == "Q" or command["desc"] == "R":
            execute = command["in"]%command["user"]
        elif command["desc"] == "qstat":
            execute = command["in"]

        # Run ls.
        child.sendline (execute)
        child.expect (COMMAND_PROMPT)
        out = child.before
        command["out"] = out

    # Now exit the remote host.
    child.sendline ('exit')
    index = child.expect([pexpect.EOF, "(?i)there are stopped jobs"])
    if index==1:
        child.sendline("exit")
        child.expect(EOF)

    return  (commands)

# parsea la salida de qstat
def qparse_outs (outs):
    for o in outs:
        lines = o["out"].split('\n')
        if o["desc"] == "Q" or o["desc"] == "R":
            jobs = []
            for l in (lines[1:len(lines) - 1]):
                w = l.split('.')
                jobs.append({"JobID":w[0]})
            o["out"] = jobs
        elif o["desc"] == "qstat":
            list_info = []
            lp = lines[1].split()
            lp[0] = lp[0] + lp[1]
            lp[4] = lp[4] + lp[5]
            lp.pop(5)
            lp.pop(1)
            list_info.append(lp)
            for l in (lines[3:len(lines) - 1]):
                lp =l.split()
                lp[0] = lp[0].split(".")[0]
                list_info.append(lp)
            info = [] 
            columns = list_info[0]

            for li in list_info[1:]:
                a = {}
                for i,e in enumerate(li):
                    a[columns[i]] = e
                info.append(a)    

            o["out"] = info 



    return outs        

if __name__ == "__main__":
    # genera el run.qsub para sheldon
    cant_nodos = 8
    qs = Qsub("UZrEquilibration", cant_nodos, "infile_equil")

    run_path = os.path.dirname(os.path.realpath(__file__))
    qs.set_run_path (run_path)
    qs.gen_file()


