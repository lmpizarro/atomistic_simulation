Modelo de Ising

crear una carpeta ./build
> cd ./build
> cmake ..
> make


=========================================================================
  PROGRAMA ISING
=========================================================================

ENTRADAS:

    1) Archivo <parametros.dat> en la misma carpeta

        El formato del archivo debe ser:

         N_R  M_R T J K_tot K_grab

        donde: 
          N_R, M_R = Tamaño de la red de spines (NxM)
          T        = Temperatura 
          J        = Constante de acoplamiento spines
          K_tot    = Cantidad de pasos Metrópolis
          K_grab   = Cada cuántos pasos se graba la 
                     energía y la magnetización
    2) Línea de comando

       ising [flag] [valor]

       donde [flag] puede ser:
      
         -t  = Temperatura
         -k  = Valor de K_tot
         -m  = Valor de K_grab
 
       Si se especifica 
       a) -m 0
       no se guardará ningún valor de energía ni de magnetización temporal
       b) -m 1
       se guardan todos los valors de energía y magnetización en cada
       paso Metrópolis (K_tot puntos por archivo)
       
       En el caso de especificar parámetros por línea de comandos, 
       éstos últimos son los que toma el programa.
       
       Ejemplo:   ./ising -t 1.5 -k 1000000 -m 10
                  
                  Corre a una temperatura de 1.5 con 1e6 pasos y graba cada 10
                  pasos

    Si existe un archivo <ultimo_estado.dat>, el programa lo lee y lo toma
    como estado inicial de la corrida

SALIDAS:
   
     Al terminar la ejecución del programa se tendrán los siguientes archivos:

     a) <val_medios.dat>
        Archivo con información sobre el valor medio y varianza de la energúa y la
        magnetización tomados sobre los K_tot pasos de la corrida
     b) <energia.dat> , <magneti.dat>
        En el caso de no haber puesto "K_grab=0" se crean estos dos archivos con
        los valores de energía y magnetización. Cada uno tendrá mod(K_tot,K_grab)+1 puntos
     c) <estado.dat>
         Configuración del último estado en formato fortran ( write(*,*) A)
         Un vector con [col_1 col_2 .... col_M]
     d) <estado_mat.dat>
         Configuración del último estado en formato matricial (NxM) valores
     e) <acpetaciones.dat>
         Archivo con la cantidad de estados aceptados en el algoritmo de Metrópolis

        
