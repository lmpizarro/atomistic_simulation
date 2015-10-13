program main_mc_pote 

    use io_parametros, only: read_parameters
    use dinmods
    use mc
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()
    call inicializacion()
    call inicia_posicion_cs()

    call metropolis()

    call finalizacion()

end program main_mc_pote
