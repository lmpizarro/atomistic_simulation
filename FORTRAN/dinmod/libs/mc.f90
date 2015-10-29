module mc
    use ziggurat
    use usozig
    use types, only: dp
    use constants
    use potenciales, only: Lenard_Jones

    use datos_problema, only : Parametros

    implicit none

    type Monte_Carlo
      real(dp),  dimension(:), allocatable    :: abs_ener
      real(dp),  dimension(:,:), allocatable  :: R
      real(dp),  dimension(:), allocatable:: obs_ener
      real(dp) :: Vol, Rho, beta
      type(Parametros) :: params
      type(Lenard_Jones) :: potencial

      integer :: r_aceptacion, kmed
    contains
      procedure :: init => inicializa
      procedure :: clear => finaliza
      procedure :: run_metropolis => metropolis_oo
      procedure :: cpc => cpc_vec
      procedure :: set_potencial => set_potencial
      procedure :: out_energ => write_energ
    end type Monte_Carlo

    private        
    public :: Monte_Carlo 

contains

  subroutine cpc_vec(this)
    class (Monte_Carlo) :: this

    this % R = this % R - this % params % gL*floor(this % R/this % params % gL)

  end subroutine

  subroutine set_potencial(this, pot)
    class (Monte_Carlo) :: this
    type(Lenard_Jones) :: pot

    this % potencial =  pot

  end subroutine


  subroutine inicializa (this, pars)
    class (Monte_Carlo) :: this
    type(Parametros) :: pars
    integer :: i, j

    this % params = pars

    this % kmed = int(this % params % gNtime / abs (this % params % gNmed)) + 1

    allocate(this % R(3,this % params % gNpart))
    allocate(this % abs_ener( this % params % gNtime))
    allocate(this % obs_ener (this % kmed))

    ! Calcula volumen y densidad
    this % Vol = this % params % gL**3
    this % Rho = this % params % gNpart / this % Vol

    !beta = 1.0 / (this % params % gT * K_BOLTZMANN)
    this % beta = 1.0 / (this % params % gT * kb)

    this % r_aceptacion = 0


    write(*,'(a)') ''
    write(*,'(a)')      '*********** PARAMETROS DERIVADOS ************'
    write(*,'(a,F8.3)') '************ Volumen               = ' , this % Vol 
    write(*,'(a,F8.4)') '************ Densidad              = ' , this % Rho
    write(*,'(a)')      '*********************************************'

    ! Inicial generador de número aleatorios
    call inic_zig()

   ! INICIALIZA Posicion aleatoria dentro de la caja 
   do i = 1, this % params % gNpart
     do j= 1, 3
        this % R(j, i) = uni() * this % params % gL 
     end do
   end do     

  end subroutine inicializa

  subroutine finaliza (this)
    class (Monte_Carlo) :: this

    ! finaliza el generador de randoms
    call fin_zig()

     ! Libera memoria
    deallocate(this % R)

  end subroutine finaliza

  !================================================================================
  ! para escribir valores de energía intermedia
  !================================================================================
  subroutine write_energ (this)
    class (Monte_Carlo) :: this
    character (len=32) :: file_energy
    logical :: es

    write(file_energy, 1000)"energy", this % params % gT, ".txt"
    print *, file_energy

    inquire(file=file_energy,exist=es)
    if ( .not. es) then
      open(unit=10,file=file_energy,status='new')
      write(10,800) this % abs_ener(10000:)
      close(10)
    end if


    800 format (E15.5)
    1000 format (A, F0.3, A)
  end subroutine write_energ

  subroutine metropolis_oo (this)
    class (Monte_Carlo) :: this
    type(Lenard_Jones):: l_j

    real(dp) :: rx, ry, rz
    real(dp) :: n_energy, deltaE
    real(dp):: pr
    integer :: iPart, i
    real(dp), dimension(3)  :: ri_vec


    n_energy = 0.0_dp
    n_energy = this % potencial  % potencial(this % R) 

    write(*,'(a)') ''
    write(*,'(a,F8.4)') '************ Energia              = ' ,n_energy 
    write(*,'(a)')      '*********************************************'

    do i=1, this % params % gNtime
      iPart = rand_int(this % params % gNpart)
      ! guardo el valor original de la energia
      rx = this % R(1, iPart) 
      ry = this % R(2, iPart)
      rz = this % R(3, iPart)

      ! actualizo el arreglo de posiciones con una nueva posicion
      this % R(1, iPart) = this % R(1, iPart) + (uni() - 0.5)
      this % R(2, iPart) = this % R(2, iPart) + (uni() - 0.5)
      this % R(3, iPart) = this % R(3, iPart) + (uni() - 0.5)
     
      ! llamo a condiciones períodicas de contorno
      call this % cpc()

      ! calculo de la variacion de energia
      !deltaE = delta_poten_lj(iPart, rx, ry, rz)
      ri_vec = (/rx, ry, rz/)
      deltaE = this % potencial  % delta_potencial(iPart, ri_vec, this % R)

      if (deltaE .lt. 0) then
          n_energy = n_energy + deltaE 
          this % r_aceptacion = this % r_aceptacion + 1
      else
          pr = exp(-this % beta * deltaE)
          if (uni() .lt. pr) then
              n_energy = n_energy + deltaE    
              this % r_aceptacion = this % r_aceptacion + 1
          else
              this % R(1, iPart) = rx
              this % R(2, iPart) = ry
              this % R(3, iPart) = rz
          endif        
      endif        
      this % abs_ener(i) = n_energy
      !if (mod(i, this % params % gNmed) .eq. 0) then
        !print *, n_energy 
      !endif

    enddo

    call this % out_energ()
    print *, this % r_aceptacion 

  endsubroutine metropolis_oo

!  subroutine cpc(l)
  
!    integer, intent(in) :: l

!    if (gR(1,l) .lt. 0) then
!      gR(1,l) = gR(1,l) + gL
!    endif        

!    if (gR(2,l) .lt. 0) then
!     gR(2,l) = gR(2,l) + gL
!    endif        

!    if (gR(3,l) .lt. 0) then
!     gR(3,l) = gR(3,l) + gL
!    endif

!    if (gR(1,l) .gt. gL) then
!     gR(1,l) = gR(1,l) - gL
!    endif        

!    if (gR(2,l) .gt. gL) then
!     gR(2,l) = gR(2,l) - gL
!    endif        

!    if (gR(3,l) .gt. gL) then
!     gR(3,l) = gR(3,l) - gL
!    endif

!  endsubroutine cpc

endmodule mc
