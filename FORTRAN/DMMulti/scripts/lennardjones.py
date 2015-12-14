import numpy as np
from numba import  double, jit, int32, void

#http://docs.continuum.io/numbapro/quickstart
@jit((double[:,:], double[:], double[:,:], double[:,:], double[:], double, double))
def lennardjones(U, box, sigma, epsilon, ind_mat, rc2, ecut):
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
            # con el indice de materiales obtengo el sgma
            # y epsilon que tengo que aplicar
            sgma = sigma [ind_mat[i], ind_mat[j]]
            epsl = epsilon [ind_mat[i], ind_mat[j]]
            print (("%f %f")%(sigma, epsilon))

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
    

@jit((double[:,:], double[:]))
def cpc_vec(U, box):
    # Periodic boundary condition
    U[0]  -= box[0] * np.rint(U[0]/box[0])
    U[1]  -= box[1] * np.rint(U[1]/box[1])
    U[2]  -= box[2] * np.rint(U[2]/box[2])


#@jit(double[:,:](double[:,:],int32, double, double, double, double[:]))
@jit
def integracion_min (U, box, rc2, ecut, gM, gNtime, gDt):

    Epot, gF, Vir = lennardjones.lennardjones(U,box, rc2, ecut)
    for i in range(gNtime):
        U = U + 0.5 * gF * gDt**2 / gM
        cpc_vec(U, box)
        Epot, gF, Vir = lennardjones.lennardjones(U,box, rc2, ecut)

    return U, Epot, gF, Vir
