module dinmods 
    use types, only: dp
    use globales, only: gT, gDt, gL, gNpart, gNtime, gR, gF, gV


    implicit none

    private
    public :: inicializacion
contains
    !===============================================================================
    ! INICIALIZA PARAMETROS
    !===============================================================================
    subroutine inicializacion()
        allocate(gR(3,gNpart))
        allocate(gV(3,gNpart))
        allocate(gF(3,gNpart))

    end subroutine inicializacion 
end module dinmods 
