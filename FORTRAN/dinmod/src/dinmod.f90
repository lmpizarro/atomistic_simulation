program main_dimod 

    use io_parametros, only: read_parameters
    use globales, only: gT, gDt, gL, gNpart, gNtime
    use constants, only: kb
    use dinmods, only: inicializacion
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()


end program main_dimod
