
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
export FC=/opt/intel/composer_xe_2013.5.192/bin/intel64/ifort
```

Hay un archivo **instalar.sh** que facilita la compilación con los flags requeridos:

```bash
mkdir build
cd build
#cmake -Doptimize=ON ..
cmake -Ddebug=OFF -Dopenmp=ON -Dverbose=OFF -Dprofile=OFF -Doptimize=ON ..
make
```

Para ver qué hace cada opción se puede mirar el archivo **CMakeList.txt** en la carpeta raíz.

Al ejecutar el script *instalar.sh* se creará una carpeta *build* donde se instalará todo el programa. Dentro de dicha carpeta se creará la carpeta *bin* donde estarán los ejecutables y scripts utilizados para postprocesar los datos.

Para ejecutar el programa, una vez dentro de la carpeta **./build/bin/** se debe hacer **./dinmod**

### Opciones del preprocesador

El programa **dinmod** tiene opciones de preprocesamiento para indicar qué cosa se calcularán durante su ejecución. Estas variables se deben especificar antes de la compilación. El archivo que se debe modificar es el **./src/control.h**

Las variables que se pueden setear son:

```
1. THERM 
  * =0 sin termostato para el ensambe microcanónico (NVE)
  * =1 termostato de Langevin para el ensamble canónico (NVT)
2. CONTROL_TEMP
  * Graba un archivo con las velocidades para luego hacer estadística y comprobar el buen funcionamiento del termostato de Langevin. **No fue probada en esta versión del programa**
3. GRABA_TRAYECTORIA
  * Graba las trayectorias de las partículas en formato VTF para luego ser visualizadas en el programa VMD
4. CORR_PAR
  * Graba un archivo con las funciones g(r) 
5. MODOS_VIB
  * Para análisis de los modos vibracionales
``` 
<! --- 6. MODOS_VIB_EQUIVALENTES
  * Para análisis de los modos vibracionales --->

B) Archivos de entrada
--------------------

El programa **dinmod** lee todos los datos de entrada del archivo **input.par**. El archivo consiste en seis líneas con números ubicados al comienzo del archivo. El resto de las líneas no son leidas por el programa y sirven para guardan configuraciónes de interés o comentarios.

El formato del archivo de entrada es el siguiente:

```Fortran
"Cantidad de especies" "Estado del sistema"
"epsilon(1)" "sigma(1)" "masa_1"
"epsilon(2)" "sigma(2)" "masa_2"
"Lado_caja" "periodo_red" "estructura_cristalina" 
"temperatura" "dt" "#_pasos" "N_med" "gamma_langevin"
"relación_especies"
```

A continuación se describen los significados de estos nombres (algunos de los datos de entrada no pueden variarse ya que hay muchas cosas aún no implementadas).

1. _Cantidad de especies_: cantidad de especies distintas presentes en el sistema. **DEBE SER 2**
2. _Estado del sistema_ : indica si se trabajará con un sólido o un líquido. **DEBE SER 1** 
3. _Epsilon(j)_: Parámetro de Lenard-Jones de la especie j
4. _Sigma(j)_: Parámetro de Lenard-Jones de la especie j
5. _Masa(j)_: Masa de la especie j
6. _Lado caja_: Lado de la caja
7. _Periodo red_: Cantidad de veces que se repite la celda primitiva en cada dimensión
8. _Estructura cristalina_: Tipo de estrictura. **DEBE SER 2** (para  FCC)
9. _Temperatura_: Temperatura del sistema
10. _dt_: Paso de integración 
11. _# pasos_: Cantidad de pasos temporales
12. _N med_: Cada cuántos puntos se realizan mediciones
13. _gamma langevin_: Parámetro para el termostato de Langevin
14. _relacion especies_:  Relación entre las distintas especies. **DEBE SER 0.5**

### Ejemplo de arcvhivo **input.par** 

```Fortran
2 1
1.0 1.0 1.0
0.8 0.9 0.6
5.95 4 2
0.5 0.001 20000 10 0.5
.5
```

C) Ejecutar con openmp
------------------------

Si se compiló el programa con la opción 'cmake -Dopenmp=ON' se ejecuta el programa
especificando primero la cantidad de threads que se quieren utilizar

```bash
export OMP_NUM_THREADS=4 && ./dinmod
```

D) Archivos de salida
-----------------------

**val_medios.dat** 
  Valores medios y desviaciones estándar de presión, temperatura y energías (potencial, cinética y total). Las primeras líneas tienen información de los datos de entrada.

**presion.dat** 
  Valores de la presión en función del tiempo. La primer linea contiene la cantidad de datos y el dt entre ellos.

**temperatura.dat**  
  Valores de la temperatura en función del tiempo. La primer linea contiene la cantidad de datos y el dt entre ellos.
ética
**energias.dat**
  Valores de las energías (potencial, cinética, y total) en función del tiempo. La primer linea contiene la cantidad de datos y el dt entre ellos.

> Las columnas son:    $E_{pot}$  $E_{cin}$   $E_{tot}$

**energia_pot_min.dat**
 Valores de la energía en función del tiempo cuando se trata de minimizarla. La primer linea contiene la cantidad de datos y el dt entre ellos

**estados.dat** 
  Posiciones, velocidades y tipo de partícula del último estado

> Formato:  $X$ $Y$ $Z$ $V_X$ $V_Y$ $V_Z$ $Tipo$

**gr.dat**
  (OPCIONAL) Sólo si se definió CORR_PAR en control.h del preprocesador. Da la función g(r) normalizada junto con la posición r. La primer linea contiene la cantidad de datos y un cero (irrelevante) 

> Formato:     $r$ $g_{11}(r)$ $g_{12}(r)$ $g_{22}(r)$ $g_{tot}(r)$                  
  
**trayectoria.vtf**
  (OPCIONAL) Sólo si se definió GRABA_TRAYECTORIA en control.h del preprocesador. Guarda las trayectorias en el formato para ser visualizado con el programa VMD.

**velocidades_control_T.dat** 
  (OPCIONAL) Sólo si se definió CONTROL_TEMP en control.h del preprocesador. Guarda las tres componentes de velocidades de una partícula cualquiera a lo largo del tiempo. Con este archivo se hace un histograma y calcula su desviación estandar y se controla que el termostato de Langevin esté funcionando. **No fue probado en esta versión del programa**

Links
----------

Fortran program examples 
[Ejemplos de Fortran](http://www.personal.psu.edu/jhm/f90/progref.html)

Fortran power
[Fortran](http://www.obliquity.com/computer/fortran/)

Programing python jobs from home
[Python Jobs](https://www.python.org/jobs/)


