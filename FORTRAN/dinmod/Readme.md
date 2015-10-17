Dinámica Molecular

crear una carpeta ./build
> cd ./build
> cmake ..
> make


formato de parametros.dat

T Npart L Dt Ntime sigma epsil

293 1000 100 .1 1000 4 1

                     Xe   Ar 
epsilon (kJ/mol)    1.77  0.997
sigma (angstroms)   4.10  4.0

donde

T: Temperatura
Npart: cantidad de partículas
L: Tamaño de la caja
Ntime: cantidad de pasos de tiempo
Dt: tamaño del paso de tiempo
sigma: potemcial de LJ
epsil: potemcial de LJ
       
!==================
! Correr con OpenMP
!==================

Si se usa gfortran, se debe compilar con el flag '-fopenmp'
Luego, para ejecutarlo se le pasa la cantidad de threads que
se quieran utilizar:

> export OMP_NUM_THREADS=2
> ./dinmod
