module globales
    
#include "control.h"

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
    real(dp)        :: gVol       ! Volumen del cubo
    real(dp)        :: gRho       ! Densidad numero de particulas N/V
   
    real(dp)        :: gPot       ! Energía potencial del sistema
    real(dp)        :: gKin       ! Energia cinetica del sistema
    real(dp)        :: gVir       ! Cálculo del virial para la presión 

    real(dp)        :: gGamma     ! Parámetro para el termostato de Langevin  


#ifdef CORR_PAR
    ! Variables utilizadas para calcular la g(r) 

    ! g(r) sin normalizar (por eso definida como integer)
    integer, dimension(:), allocatable   :: gCorr_par 

    integer                 :: gNhist     ! Cantidad de bines de la g(r)
    integer                 :: gNgr       ! Contador para saber cuántas veces se acumuló la g(r)
    real(dp)                :: gDbin      ! Ancho del bin de la g(r)    
   
#endif


end module globales
