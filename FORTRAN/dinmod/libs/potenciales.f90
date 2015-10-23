module potenciales 
  use types, only: dp
  !use globales, only: gF, gV, gR
  use datos_problema, only : Parametros

  implicit none


  type Lenard_Jones
    real(dp)        :: rc2       ! Cuadrado del radio de corte
    real(dp)        :: pot_cut   ! Potencial L-J en el radio de corte
    type(Parametros) :: params
 
  contains
    procedure :: calc_radio_corte => corta_desplaza_pote
    procedure :: delta_potencial => delta_poten_lj_vec
    procedure :: potencial => poten_lj_vec
    procedure :: kernel => kernel_pot
    procedure :: init => set_params
  end type Lenard_Jones

  private
  public :: Lenard_Jones 

contains

  subroutine set_params(this, pars)   
    class (Lenard_Jones) :: this
    type(Parametros) :: pars

    this % params = pars

    write(*,'(a)') ''
    write(*,'(a)')      '********  Parámetros en el Potencial ********'
    write(*,'(a,F8.4)') '************ Sigma              = ' , this % params % gSigma 
    write(*,'(a,F8.4)') '************ Sigma              = ' , this % params % gEpsil 
    write(*,'(a)')      '*********************************************'


  end subroutine set_params   


  !===============================================================================
  ! RADIO DE CORTE 
  !===============================================================================
  ! Define el radio de corte y calcula el desplazamiento del potencial de L-J 
  subroutine corta_desplaza_pote(this)   
    class (Lenard_Jones) :: this

    ! Radio de corte fijo
    this % rc2 = (2.5_dp)**2                
    ! Potencial de L-J evaluado en el radio de corte
    this % pot_cut = 4.0_dp*this % params % gEpsil*( 1.0_dp/(this % rc2**6) - 1.0_dp/(this % rc2**3) ) 
    
  end subroutine corta_desplaza_pote

  !
  ! Calculo del potencial en forma vectorial
  ! 
  function poten_lj_vec (this, R) result(Pot)
    class (Lenard_Jones) :: this
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp)                :: Pot
    integer                 :: i,j       
    real(dp),  dimension(:,:) :: R

    call this % calc_radio_corte()

    ! Paso a trabajar distancias en unidades de sigma
    R = R/this % params % gSigma
    this % params % gL = this % params % gL/this % params % gSigma

    Pot = 0.0_dp

    do i = 1, this % params % gNpart - 1       
      do j = i+1, this % params % gNpart
        rij_vec = R(:,i) - R(:,j)               ! Distancia vectorial
        ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
        ! Siempre en distancias relativas de sigma
        rij_vec = rij_vec - this % params % gL*nint(rij_vec/this % params % gL)
        r2ij   = dot_product( rij_vec , rij_vec )    ! Cuadrado de la distancia
        if ( r2ij < this % rc2 ) then               
          r2in = 1.0_dp/r2ij                         ! Inversa al cuadrado
          r6in = r2in**3                             ! Inversa a la sexta
          Pot     = Pot + r6in * ( r6in - 1.0_dp)    ! Energía potencial
        end if
      enddo
    enddo

    ! Constantes que faltaban en el potencial
    ! Agrego el desplazamiento del potencial considerando la cantidad de
    ! pares con que se obtuvo la energía potencial N(N-1)/2
    Pot =  4.0_dp * this % params % gEpsil * Pot - this % params % gNpart*(this % params % gNpart-1)* this % pot_cut/2.0_dp  

    ! Se vuelven a pasar a las coordenadas absolutas
    R = R*this % params % gSigma
    this % params % gL = this % params % gL*this % params % gSigma

  endfunction poten_lj_vec

  !
  !  calculo del delta de potencial en forma vectorial
  !
  function delta_poten_lj_vec  (this, i, ri_vec, R) result(deltaPot)
    class (Lenard_Jones) :: this
    real(dp), intent(in), dimension(:,:)  :: R

    ! la particula con la que se hace le calculo
    integer, intent(in) :: i    
    ! la posición original de la partícula 
    real(dp), intent(in), dimension(3)  :: ri_vec   
    integer :: j
    real(dp) :: deltaPot
    real(dp) :: vo, vn

    deltaPot = 0.0_dp

    vo = 0.0_dp
    vn = 0.0_dp
    do j = 1, i - 1
      vo     = vo + this % kernel(j, ri_vec, R)   ! Energía potencial
      vn     = vn + this % kernel(j, R(:,i), R)   ! Energía potencial
    enddo

    do j = i + 1, this % params % gNpart 
      vo     = vo + this % kernel(j, ri_vec, R)   ! Energía potencial
      vn     = vn + this % kernel(j, R(:,i), R)   ! Energía potencial
    enddo

    deltaPot =  4.0_dp * this % params % gEpsil * (vn - vo)   

  endfunction delta_poten_lj_vec

  !
  ! una funcion para simplificar los calculos del delta de potencial
  !
  function kernel_pot (this, j, ri_vec, R) result (Pot)
    class (Lenard_Jones) :: this
    real(dp), intent(in), dimension(3)  :: ri_vec
    integer, intent(in) :: j
    real(dp)                :: r2ij      ! Módulo cuadrado de la distancia rij
    real(dp)                :: r2in,r6in ! Inversa distancia rij a la 2 y 6
    real(dp)                :: Pot
    real(dp), dimension(3)  :: rij_vec   ! Distancia vectorial entre i y j
    real(dp),  dimension(:,:) :: R

    Pot = 0

    rij_vec = ri_vec - R(:,j)               ! Distancia vectorial
    ! Si las partícula está a más de gL/2, la traslado a r' = r +/- L
    ! Siempre en distancias relativas de sigma
    rij_vec = rij_vec - this % params % gL*nint(rij_vec/this % params % gL)
    r2ij   = dot_product( rij_vec , rij_vec )    ! Cuadrado de la distancia
    if ( r2ij < this % rc2 ) then               
       r2in = 1.0_dp/r2ij                         ! Inversa al cuadrado
       r6in = r2in**3                             ! Inversa a la sexta
       Pot = r6in * ( r6in - 1.0_dp)    ! Energía potencial
    end if
  endfunction kernel_pot

end module potenciales 
