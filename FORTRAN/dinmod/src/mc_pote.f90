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

    ! Lee los datos del problema
    call params % leer()
    ! inicializa montecarlo 
    call mc % init(params)
    ! inicializa el potencial
    call lj % init(params)

    ! setea el potencial en Monte Carlo
    call mc % set_potencial(lj)

    ! Corre metropolis
    call mc % run_metropolis()

    call mc % out_energ()
    call mc % out_presion()
    print *, mc % r_aceptacion 



    call mc % clear()

end program main_mc_pote
