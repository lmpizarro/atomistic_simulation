
Dinámica Molecular para un sistema binaro
=========================================

El programa resuelve un sistema de dos partículas que interactúan con un potencial de L-J.
Está pensado para trabajar con sistemas sólidos, pero puede adaptarse fácilmente para otra
situación.

A) Instalación
-----------

Para compilar el programa se requiere tener instalado el paquete *cmake*. En el archivo CMakeList están todos los comandos que serán pasados al compilador. Se puede elegir el compilador a utilizar (gfortran, ifrot, etc) definiendo la variable de entorno *FC*:
```bash
# por ejemplo en sheldon se usa 
# export FC=/opt/intel/composer_xe_2013.5.192/bin/intel64/ifort
```

Hay un archivo *instalar.sh* que facilita la compilación con los flags requeridos:

```bash
mkdir build
cd build
#cmake -Doptimize=ON ..
cmake -Ddebug=OFF -Dopenmp=ON -Dverbose=OFF -Dprofile=OFF -Doptimize=ON ..
make
```

Para ver qué hace cada opción se puede mirar el archivo *CMakeList.txt* en la carpeta raíz.

Al ejecutar el script *instalar.sh* se creará una carpeta *build* donde se instalará todo el programa. Dentro de dicha carpeta se creará la carpeta *bin* donde estarán los ejecutables y scripts utilizados para postprocesar los datos.

Para ejecutar el programa, una vez dentro de la carpeta *./build/bin/* se debe hacer *./dinmod*

### Opciones del preprocesador

El programa *dinmod* tiene opciones de preprocesamiento para indicar qué cosa se calcularán durante su ejecución. Estas variables se deben especificar antes de la compilación. El archivo que se debe modificar es el *./src/control.h

Las variables que se pueden setear son:

```c
1. THERM 
  *=0 sin termostato para el ensambe microcanónico (NVE)
  * =1 termostato de Langevin para el ensamble canónico (NVT)
2. CONTROL_TEMP
  * Graba un archivo con las velocidades para luego hacer estadística y comprobar el buen funcionamiento del termostato de Langevin
3. GRABA_TRAYECTORIA
  * Graba las trayectorias de las partículas en formato *vtf* para luego ser visualizadas en el programa *VMD*
4. CORR_PAR
  * Graba un archivo con las funciones g(r) 
5. MODOS_VIB
  * Para análisis de los modos vibracionales
6. MODOS_VIB_EQUIVALENTES
  * Para análisis de los modos vibracionales
```

1. THERM 
  *=0 sin termostato para el ensambe microcanónico (NVE)
  * =1 termostato de Langevin para el ensamble canónico (NVT)
2. CONTROL_TEMP
  * Graba un archivo con las velocidades para luego hacer estadística y comprobar el buen funcionamiento del termostato de Langevin
3. GRABA_TRAYECTORIA
  * Graba las trayectorias de las partículas en formato *vtf* para luego ser visualizadas en el programa *VMD*
4. CORR_PAR
  * Graba un archivo con las funciones g(r) 
5. MODOS_VIB
  * Para análisis de los modos vibracionales
6. MODOS_VIB_EQUIVALENTES
  * Para análisis de los modos vibracionales

B) Archivos de entrada
--------------------

El programa *dinmod* lee todos los datos de entrada del archivo *input.par*. El archivo consiste en seis líneas con números ubicados al comienzo del archivo. El resto de las líneas no son leidas por el programa y sirven para guardan configuraciónes de interés o comentarios.

El formato del archivo de entrada es el siguiente:

_<Cantidad de especies>_ _<Estado del sistema>_
_epsilon(1)_ _sigma(1)_ _masa_1_
_epsilon(2)_ _sigma(2)_ _masa_2_
_temperatura_ _dt_ _# pasos_ _N_med_ _gamma_langevin_
_relación_especies_

### Ejemplo de programa

```Fortran
2 1
1.0 1.0 1.0
0.8 0.9 0.6
5.95 4 2
0.5 0.001 20000 10 0.5
.5
```

Formato del archivo input.par

cant_especies liqsol
epsilon_1 sigma_1 masa_1
epsilon_2 sigma_2 masa_2
.
.
.
epsilon_n sigma_n masa_n
Si es gas -> Lcaja densidad temperatura delta_t  total pasos pasos_medicion gamma
Si es solido -> Lcaja cantidad_celdas tipo_cubica temperatura delta_t  total pasos pasos_medicion gamma
pc1
pc2
.
.
.
pc_n-1

4
1.1E-2 4.1 1.1
1.2E-2 4.2 1.2
1.3E-2 4.3 1.3
1.4E-2 4.4 1.4
5 .3  1.2 0.001 10000 100 0.5
.2
.3
.4

liqsol 0   posiciones aleatorias
liqsol 1   posiciones cristalinas

cantidad_celdas en x, y, z
tipo_cubica 0: simple, 1: centrada en cuerpo, 2: centrada en caras

<!-- fortran program examples -->
Ejemplos [Fortran](http://www.personal.psu.edu/jhm/f90/progref.html)

fortran power
Referencia [Fortran](http://www.obliquity.com/computer/fortran/)

programing python jpbs from home
[Python Jobs](https://www.python.org/jobs/)

---------------------------------------------------------------------
- EJECUTAR CON OPENMP
---------------------------------------------------------------------

Si se compiló el programa con la opción 'cmake -Dopenmp=ON' se ejecuta el programa
especificando primero la cantidad de threads que se quieren utilizar

 > export OMP_NUM_THREADS=4 && ./dinmod

