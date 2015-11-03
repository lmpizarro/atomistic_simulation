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
4) Cálculo de la función g(r)
   #define/undef CORR_PAR -> activa opciones de cálculo de la función
        de correlación de pares g(r). El cálculo principal se hace
        dentro del loop de fuerzas (paralelizado  y en serie). También
        hay porciones de codigo en la inicializacion y finalización
        del programa. Se basa en el algoritmo de la pg. 86 de Frenkesl
        El resultado se guara en el archivo <gr.dat> con el formato:      
                     r(1)  gr(1)
                     r(2)  gr(2)
                      ....... 
                     r(N)  gr(N)
        El valor de N (cantidad de histogramas) está hardcodeado en la 
        inicialización del programa.
            
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

----------------------------------------------------------------------
- ARCHIVOS DE SALIDA
----------------------------------------------------------------------

 <val_medios.dat>  : Valores medios y desviaciones estándar de presión, temperatura 
                     y energías (potencial, cinética y total).
                     Las primeras líneas tienen información de los datos de entrada

 <presion.dat>     : Valores de la presión en función del tiempo
                     La primer linea contiene la cantidad de datos y el dt entre ellos

 <temperatura.dat> : Valores de la temperatura en función del tiempo
                     La primer linea contiene la cantidad de datos y el dt entre ellos

 <energias.dat>    : Valores de las energías (pot, cin, y total) en función del tiempo
                     La primer linea contiene la cantidad de datos y el dt entre ellos
                     Las columnas son:    E_POT   E_CIN   E_TOT

 <energia_pot_min.dat> : Valores de la energía en función del tiempo cuando se trata
                         de minimizarla

 <estados.dat>     : Posiciones y velocidades del último estado
                     Formato:   X Y Z V_X V_Y V_Z

 <gr.dat>          : (OPCIONAL) Sólo si se definió CORR_PAR en control.h del preprocesador
                     Da la función g(r) normalizada junto con la posición r.
                     La primer linea contiene la cantidad de datos y un cero (irrelevante) 
                     Formato:     r  g(r)                   
  
 <trayectoria.vtf>  : (OPCIONAL) Sólo si se definió GRABA_TRAYECTORIA en control.h del
                      preprocesador. Guarda las trayectorias en el formato para ser 
                      visualizado con el programa VMD.

<velocidades_control_T.dat> : (OPCIONAL) Sólo si se definió CONTROL_TEMP en control.h del
                              preprocesador. Guarda las tres componentes de velocidades
                              de una partícula cualquiera a lo largo del tiempo. Con este
                              archivo se hace un histograma y calcula su desviación estandar
                              y se controla que el termostato de Langevin esté funcionando.
