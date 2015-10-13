program main_dimod 

    use io_parametros, only: read_parameters
    use globales, only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM
    use constants, only: kb
    use dinmods, only: inicializacion, inicia_posicion_cs, finalizacion, fuerza
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()
    call inicializacion()

    call inicia_posicion_cs ()

    call fuerza()

    call finalizacion()

end program main_dimod
