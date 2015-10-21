module globales
    
    use types, only: dp

    implicit none

    ! gT: Temperatura de la corrida
    ! gDt: Paso de tiempo de la corrida
    ! gL: longitud de un lado del cubo
    real(dp) :: gT, gDt, gL
    ! gNpart: cantidad de partículas del sistema
    ! gNtime: cantidad de pasos de  tiempo
    ! gNmed : cantidad de pasos entre mediciones
    integer :: gNpart, gNtime, gNmed
    ! sigma: 
    ! epsil:
    ! parametros del potencial de LJ
    real(dp) :: gSigma, gEpsil
    ! gR:  posicion de la particula
    ! gF:  fuerza entre particulas particula
    ! gV:  velocidad de la particula
    ! gM:  masa de la partícula
    real(dp),  dimension(:,:), allocatable  :: gR, gF, gV
    real(dp)                                :: gM
    real(dp),  dimension(:), allocatable    :: gAbs_ener

    real(dp)        :: gRc2       ! Cuadrado del radio de corte
    real(dp)        :: gPot_cut   ! Potencial L-J en el radio de corte
   
end module globales
