module globales
    use types, only: dp

    implicit none

    ! gT: Temperatura de la corrida
    ! gDt: Paso de tiempo de la corrida
    ! gL: longitud de un lado del cubo
    real(dp) :: gT, gDt, gL
    ! gNpart: cantidad de partículas del sistema
    ! gNtime: cantidad de pasos de  tiempo
    integer :: gNpart, gNtime
    ! sigma: 
    ! epsil:
    ! parametros del potencial de LJ
    real(dp) :: sigma, epsil
    ! gR:  posicion de la particula
    ! gF:  fuerza entre particulas particula
    ! gV:  velocidad de la particula
    real(dp),  dimension(:,:), allocatable:: gR, gF, gV
    real(dp),  dimension(:), allocatable:: obs_ener

end module globales
