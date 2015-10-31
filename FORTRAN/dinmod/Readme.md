Dinámica Molecular

crear una carpeta ./build
> cd ./build
> cmake ..
> make

 Para pasarle opciones de compilación al cmake:
 cmake -D[opc]=ON ..
 donde [opc] pueden ser:
   debug
   profile
   optimize
   openmp
   verbose

  Por default, las opciones del falg son -cpp -std=f2003 -fbacktrace
  Todo esto se puede modificar desde el archivo CMakelists.txt

  Ejemplo:
         >    mkdir build && cd build
         >    cmake -Ddebug=ON ..
         >    make

 Los flags se van acumulando si se ponen varias opciones simultáneas

 Sin OpenMC hacer:

         >     cmake -Dopenmp=OFF ..

--------------------------------------------------------------------
 Con el preprocesador se controlas los siguientes parámetros
 definidos en el archivo <control.h>: 

1) Termostato
   #define THERM 0 -> No se utiliza termostato ( E constante)
   #define THERM 1 -> Termostato de langevin
2) Corroborar la temperatura
   #define/undef CONTROL_TEMP -> Activa opciones para guardar las
        componentes de temperatura de una partícula arbitraria. Con 
        ellas se puede obtener un histograma y verificar que
        std(v_i) = sqrt(T)
        * Se escribe al archivo <velocidades_control_T.dat>
3) Guardar trayectorias de las partículas para visualizar con VMD
   #define/undef GRABA_TRAYECTORIA -> Activa opciones para guardar
        las posiciones de las partículas con el formato para ser
        visualizados con el programa VMD.
        * Se escribe al archivo <trayectoria.vtf>

---------------------------------------------------------------------

formato de parametros.dat

T Npart L Dt Ntime sigma epsil Masa Nmed
Gamma (para el termostato de Langevin)
293 1000 100 .1 1000 4 1 1 100
0.5

                     Xe   Ar 
epsilon (kJ/mol)    1.77  0.997 (1.83592*10-2 1.03413*10-2 en eV/at)
sigma (angstroms)   4.10  4.0

donde

<<<<<<< HEAD
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
=======
T      : Temperatura
Npart  : cantidad de partículas
L      : Tamaño de la caja
Ntime  : cantidad de pasos de tiempo
Dt     : tamaño del paso de tiempo
sigma  : potemcial de LJ
epsil  : potemcial de LJ
M      : Masa de la partícula
Nmed   : abs(Nmed) Cada cuántos pasos se realizan mediciones       
         Se Nmed > 0 se graban datos temporales (energía,persión temperatura)
         Si Nmed <= 0 no se graba ningún datos temporales (sí los
         promedios y desviaciones estandar)

---------------------------------------------------------------------
- EJECUTAR CON OPENMP
---------------------------------------------------------------------

Si se compiló el programa con la opción 'cmake -Dopenmp=ON' se ejecuta el programa
especificando primero la cantidad de threads que se quieren utilizar

 > export OMP_NUM_THREADS=4 && ./dinmod

>>>>>>> origin/master
