program main_mc_pote 

    use io_parametros, only: read_parameters
    use dinmods
  

    implicit none

    ! Lee los datos necesario
    call read_parameters()
    call inicializacion()


    call finalizacion()

end program main_mc_pote
