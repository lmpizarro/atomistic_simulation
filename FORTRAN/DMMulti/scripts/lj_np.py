
# coding: utf-8

# In[106]:

# http://combichem.blogspot.com.ar/2013/04/fun-with-numba-numpy-and-f2py.html

import numpy as np

Na = 4
Nb = 5

L = 50
import numpy
from numba import double, jit

rc2 = 6.25
rc2i=1.0/rc2
rc6i=rc2i*rc2i*rc2i
ecut=rc6i*(rc6i-1.0)


#@jit(argtypes=[double[:,:], double[:]])
def lennardjones(U, box):

    # Can't use (ndim, npart) = numpy.shape(U)
    # with Numba. No unpacking of tuples.

    ndim = len(U)
    npart = len(U[0])

    F = numpy.zeros((ndim, npart))

    Epot = 0.0
    Vir = 0.0

    for i in range(npart):
        for j in range(npart):
            if i > j:

                X  = U[0, j] - U[0, i]
                Y  = U[1, j] - U[1, i]
                Z  = U[2, j] - U[2, i]

                # Periodic boundary condition
                X  -= box[0] * numpy.rint(X/box[0])
                Y  -= box[1] * numpy.rint(Y/box[1])
                Z  -= box[2] * numpy.rint(Z/box[2])

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
    


# In[107]:

import numpy
from forces import lennardjones

# Initialize simulations (code not shown)
temperature = 2.0
n_atoms = 100
rho = 0.01

U = initialize_positions(n_atoms, rho)
box = initialize_box(n_atoms, rho)
V = initialize_velocities(U, temperature)

# Start simulation

n_steps = 10000
dt = 0.001

(epot,F,Vir) = lennardjones(U,box)

for i in range(n_steps):
    U += V * dt + 0.5 * F * dt * dt
    F0 = F[:]

    (epot, F, Vir) = lennardjones(U, box)

    V += 0.5 * (F + F0) * dt


# In[ ]:



