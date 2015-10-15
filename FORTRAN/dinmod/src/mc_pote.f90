program main_mc_pote 

    use types
    use io_parametros, only: read_parameters
    use dinmods
    use mc
    use potenciales
  

    implicit none

    real(dp) :: pot 

    ! Lee los datos necesario
    call read_parameters()
    call inicializacion()
    call inicia_posicion_rn()

    !call metropolis()

    pot = poten_lj_vec() 

    print *, "potenciales lj vec", pot

    call finalizacion()

end program main_mc_pote
