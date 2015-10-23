program main_mc_pote 

    use types
   ! use dinmods
    use mc, only : Monte_Carlo
    use datos_problema, only : Parametros
    use potenciales, only : Lenard_Jones
  

    implicit none

    type(Parametros) :: params
    type(Monte_Carlo) :: mc
    type(Lenard_Jones) :: lj

    ! Lee los datos necesario
    call params % leer()
    call mc % init(params)
    call lj % init(params)
    call mc % set_potencial(lj)

    call mc % run_metropolis()


    call mc % clear()

end program main_mc_pote
