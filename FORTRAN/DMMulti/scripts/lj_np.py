import numpy as np
from numba import double, jit

# http://combichem.blogspot.com.ar/2013/04/fun-with-numba-numpy-and-f2py.html

Na = 80
Nb = 128
Npart = Na + Nb

a=np.zeros(Na+Nb)
a[Na:Na+Nb] = 1

rc2 = 4
ecut = -1

for i in range (Npart-1):
    for j in range(i+1, Npart):
	#print i, j, a[i],a[j]
        pass

#http://docs.continuum.io/numbapro/quickstart
@jit((double[:,:], double[:]))
def lennardjones(U, box):
    # Can't use (ndim, npart) = numpy.shape(U)
    # with Numba. No unpacking of tuples.

    ndim = len(U)
    npart = len(U[0])

    F = np.zeros((ndim, npart))

    Epot = 0.0
    Vir = 0.0

    for i in range(npart - 1):
        for j in range(i+1, npart):
                
        	#print (("%d %d %d %d")%(i, j, a[i],a[j]))

            X  = U[0, j] - U[0, i]
            Y  = U[1, j] - U[1, i]
            Z  = U[2, j] - U[2, i]

                # Periodic boundary condition
            X  -= box[0] * np.rint(X/box[0])
            Y  -= box[1] * np.rint(Y/box[1])
            Z  -= box[2] * np.rint(Z/box[2])

            # Distance squared
            r2 = X*X + Y*Y + Z*Z
            if(r2 < rc2):
                r2i = 1.0 / r2
                r6i = r2i*r2i*r2i
                Epot = Epot + r6i*(r6i-1.0) - ecut

                ftmp = 48. * r6i*(r6i- 0.5) * r2i

                F[0, i] -= ftmp * X
                F[1, i] -= ftmp * Y
                F[2, i] -= ftmp * Z
                F[0, j] += ftmp * X
                F[1, j] += ftmp * Y
                F[2, j] += ftmp * Z
                Vir += ftmp

    Epot = Epot * 4.0

    return Epot, F, Vir
    


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
box = np.ones(3)
U = initialize_positions(Npart, L)
V = InitVelocities(Npart, temperature)


Epot, F, Vir = lennardjones(U,box)

print Epot, F, Vir

