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

---------------------------------------------------------------------

formato de parametros.dat

T Npart L Dt Ntime sigma epsil M Nmed

293 1000 100 .1 1000 4 1 1 100

                     Xe   Ar 
epsilon (kJ/mol)    1.77  0.997
sigma (angstroms)   4.10  4.0

donde

T      : Temperatura
Npart  : cantidad de partículas
L      : Tamaño de la caja
Ntime  : cantidad de pasos de tiempo
Dt     : tamaño del paso de tiempo
sigma  : potemcial de LJ
epsil  : potemcial de LJ
M      : Masa de la partícula
Nmed   : Cada cuántos pasos se realizan mediciones       


---------------------------------------------------------------------
- EJECUTAR CON OPENMP
---------------------------------------------------------------------

Si se compiló el programa con la opción 'cmake -Dopenmp=ON' se ejecuta el programa
especificando primero la cantidad de threads que se quieren utilizar

 > export OMP_NUM_THREADS=4 && ./dinmod

