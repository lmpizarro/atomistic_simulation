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
    call mc % init(params)

    !call mc % run_metropolis()


    call mc % clear()

end program main_mc_pote
