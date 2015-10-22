program main_mc_pote 

    use types
    use dinmods
    use mc, only : Monte_Carlo
    use datos_problema, only : Parametros
  

    implicit none

    real(dp) :: pot 
    type(Parametros) :: params
    type(Monte_Carlo) :: mc

    ! Lee los datos necesario
    call params % leer()
    call inicializacion()
    call inicia_posicion_rn()

    call mc % run_metropolis(params)

    !pot = poten_lj_vec() 

    !print *, "potenciales lj vec", pot

    !call finalizacion()

end program main_mc_pote
