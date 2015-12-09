import numpy as np
import lennardjones

# http://combichem.blogspot.com.ar/2013/04/fun-with-numba-numpy-and-f2py.html

Na = 8
Nb = 12
Npart = Na + Nb

ind_mat=np.zeros(Na+Nb)
ind_mat[Na:Na+Nb] = 1

# espilon sigma masa
el_1 = np.asarray([1, 1, 1])
el_2 = np.asarray([1.1, 1.3, 2])

#calculo de las matrices sigma epsilon
sigma12=sigma21= (el_1[1] + el_2[1] ) / 2
epsilon12=epsilon21= np.sqrt((el_1[0] * el_2[0] )) 

sigma = np.array((el_1[1], sigma12, sigma21, el_2[1]))
sigma = sigma.reshape(2,2) / el_1[1]
epsilon = np.array((el_1[0], epsilon12, epsilon21, el_2[0]))
epsilon = epsilon.reshape(2,2) / el_1[0]

# calculo de los rc y los ecut
rc = 2.5 * sigma[0,0]
rc2 = np.power(rc, 2)
el_1 = np.hstack((el_1, rc2))
ecut = 4 * (np.power(rc, -12) - np.power(rc, -6))
el_1 = np.hstack((el_1, ecut))
print "ecut", el_1 

rc = 2.5 * sigma[1,1]
rc2 = np.power(rc, 2)
el_2 = np.hstack((el_2, rc2))
ecut = 4 * (np.power(rc, -12) - np.power(rc, -6))
el_2 = np.hstack((el_2, ecut))
print "ecut", el_2 

rc2 = 4.0
ecut = -1.0
gDt = 0.001
gM = 1.0


def initialize_box(n_atoms, rho):
  # rho = natoms / l3	
  return(np.rint (np.power(n_atoms/rho, 1.0/3.0)))

def RescaleVelocities(Vel, T):
    """Rescales velocities in the system to the target temperature.
Input:
    Vel: (N,3) array of atomic velocities
    T: target temperature
Output:
    Vel: same as above
"""
    #recenter to zero net momentum (assuming all masses same)
    Vel = Vel - Vel.mean(axis=0)
    #find the total kinetic energy
    KE = 0.5 * np.sum(Vel * Vel)
    #find velocity scale factor from ratios of kinetic energy
    VScale = np.sqrt(1.5 * len(Vel) * T / KE)
    Vel = Vel * VScale
    return Vel  


def InitVelocities(N, T):
    """Returns an initial random velocity set.
Input:
    N: number of atoms
    T: target temperature
Output:
    Vel: (N,3) array of atomic velocities
"""
    Vel = np.random.rand(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel

def InitPositions(N, L):
    """Returns an array of initial positions of each atom,
placed on a cubic lattice for convenience.
Input:
    N: number of atoms
    L: box length
Output:
    Pos: (N,3) array of positions
"""
    #make the position array
    Pos = np.zeros((N,3), float)
    #compute integer grid # of locations for cubic lattice
    NLat = int(N**(1./3.) + 1.)
    #make an array of lattice sites
    r = L * (np.arange(NLat, dtype=float)/NLat)
    #loop through x, y, z positions in lattice until done
    #for every atom in the system
    i = 0
    for x in r:
        for y in r:
            for z in r:
                Pos[i] = np.array([x,y,z], float)
                i += 1
                #if done placing atoms, return
                if i >= N:
                    return Pos
    return Pos

def initialize_positions(N, L):
    """Returns an array of initial positions of each atom,
placed on a cubic lattice for convenience.
Input:
    N: number of atoms
    L: box length
Output:
    Pos: (N,3) array of positions
"""
    #make the position array
    Pos = np.random.rand(N*3).reshape(3,N) * L


    return Pos


rho = .05
temperature = 1.0

L = initialize_box(Npart, rho)
box = np.ones(3) * L
U = initialize_positions(Npart, L)
V = InitVelocities(Npart, temperature)


Epot, F, Vir = lennardjones.lennardjones(U,box, sigma, epsilon, ind_mat, rc2, ecut)

print Epot, F, Vir

#U, Epot, F, Vir = lennardjones.integracion_min(U,box, rc2, ecut, 10, gDt)

# minimiza la energia
for i in range(100):
    U = U + 0.5 * F * gDt**2 / gM
    lennardjones.cpc_vec(U, box)

print "despues"
print Epot, F, Vir

