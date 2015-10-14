program main_dimod 

    use io_parametros,    only: read_parameters
    use globales,         only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM
    use constants,        only: kb
    use dinmods,          only: inicializacion, finalizacion, &
                                fuerza, integracion_min, integracion 
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()
 
    call inicializacion()

    !call integracion_min()

    call integracion()

    call finalizacion()

end program main_dimod
