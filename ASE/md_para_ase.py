import numpy as np 
from asap3 import Atoms, EMT, units
from ase.visualize.primiplotter import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.langevin import Langevin

from ase.io.trajectory import Trajectory

from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#/home/pizarro/opt/openmpi-1.8.1-gcc/bin/mpirun -np 4 /home/pizarro/python/bin/python ./Melting.py

Temperatures = []
Tinic = 200
Tfin = 2000
Step_T = 200
for t in range(Tinic,Tfin + Step_T, Step_T):
    Temperatures.append(t)

Cant_Puntos = len(Temperatures) 
Puntos_Por_Nodo = Cant_Puntos / size + 1

# Create the atoms
atoms = FaceCenteredCubic(size=(2,2,2), symbol="Cu", pbc=False)

# Associate the EMT potential with the atoms
atoms.set_calculator(EMT())

# How many steps at each temperature
nsteps = 5000

# Make the Langevin dynamics module
dyn = Langevin(atoms, 5*units.fs, units.kB*Temperatures[0], 0.002)

# Some functions for calculating the actual temperature, energy, ...

def temperature(a):
    return 2.0/3.0 * a.get_kinetic_energy() / (len(a) * units.kB)

def etotal(a):
    return (a.get_kinetic_energy() + a.get_potential_energy()) / len(a)

def calcula_dyn (t, dyn):
    dyn.set_temperature(units.kB*t)
    data = np.empty(5, dtype=np.float64) 
    for i in range(nsteps/100):
        dyn.run(100)
        #data_step = {"E_total":  "T":  "goal": "step": "n_step": }
        data_step = np.asarray([etotal(atoms), temperature(atoms), t, i, nsteps/100], dtype=np.float64)
        data = np.vstack((data, data_step))
    data = np.delete(data, 0, 0)
    return(data)

# para acumular los resultados parciales
data =np.zeros(5, dtype=np.float64) 

#buffer para recibir los datos de los otros nodos
rec_buf = np.empty(( 151, 5), dtype=np.float64) 

# para acumular los resultados de todos los nodos 
acum_data =np.empty(5, dtype=np.float64) 

Ti = rank * Puntos_Por_Nodo
Tf = (rank + 1)  * Puntos_Por_Nodo
temps =  Temperatures[Ti: Tf]

for ts in temps:
    tmp = calcula_dyn (ts, dyn)
    data  = np.vstack((data, tmp))

print rank, temps, np.shape(data)
#print data

if rank == 0:
    acum_data = np.vstack((acum_data, data))
    for i in range(1,size):
        comm.Recv(rec_buf, ANY_SOURCE)
        print "shape: ", np.shape(rec_buf)
        acum_data = np.vstack((acum_data, rec_buf))
        pass
else:
    comm.Send(data, dest=0)

mask = np.all(np.abs(acum_data) < 1e-6, axis=0)

#acum_data = acum_data[~np.all(acum_data == 0, axis=0)]

if rank == 0:
    print "fin: ", np.shape(acum_data), np.shape(mask)
    print acum_data

