program main_dimod 

    use io_parametros,    only: read_parameters
    use globales,         only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM
    use constants,        only: kb
    use dinmods,          only: inicializacion, finalizacion, &
                                integracion_min, integracion 
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()
    ! Inicializa parámetros 
    call inicializacion()
    ! Busca el mínimo de energía
    call integracion_min()
    ! Integración de la dinámica
    call integracion()
    ! Finalización del programa
    call finalizacion()

end program main_dimod
