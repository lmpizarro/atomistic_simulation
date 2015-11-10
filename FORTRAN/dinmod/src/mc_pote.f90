program main_mc_pote 

    use types
   ! use dinmods
    use mc, only : Monte_Carlo
    use datos_problema, only : Parametros
    use potenciales, only : Lenard_Jones
  

    implicit none

    type(Parametros) :: params
    type(Monte_Carlo) :: mc_l
    type(Lenard_Jones) :: lj

    ! Lee los datos del problema
    call params % leer()
    ! inicializa montecarlo 
    call mc_l % init(params)
    ! inicializa el potencial
    call lj % init(params)

    ! setea el potencial en Monte Carlo
    call mc_l % set_potencial(lj)

    ! Corre metropolis
    call mc_l % run_metropolis()

    call mc_l % out_energ()
    call mc_l % out_presion()
    print *, mc_l % r_aceptacion 



    call mc_l % clear()

end program main_mc_pote
