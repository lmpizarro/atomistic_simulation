module mc
    use ziggurat
    use usozig
    use types, only: dp
    use constants
    use potenciales, only: Lenard_Jones
    use dinmods

    use datos_problema, only : Parametros

    implicit none

    type Monte_Carlo
      real(dp),  dimension(:), allocatable    :: abs_ener
      real(dp),  dimension(:,:), allocatable  :: R
      real(dp) :: Vol, Rho
      type(Parametros) :: params
    contains
      procedure :: init => inicializa
      procedure :: clear => finaliza
      procedure :: run_metropolis => metropolis_oo
    end type Monte_Carlo

    private        
    public :: Monte_Carlo 

contains

  subroutine inicializa (this, pars)
    class (Monte_Carlo) :: this
    type(Parametros) :: pars

    this % params = pars

    allocate(this % R(3,this % params % gNpart))


    ! Calcula volumen y densidad
    this % Vol = this % params % gL**3
    this % Rho = this % params % gNpart / this % Vol

    write(*,'(a)') ''
    write(*,'(a)')      '*********** PARAMETROS DERIVADOS ************'
    write(*,'(a,F8.3)') '************ Volumen               = ' , this % Vol 
    write(*,'(a,F8.4)') '************ Densidad              = ' , this % Rho
    write(*,'(a)')      '*********************************************'

    ! Inicial generador de número aleatorios
    call inic_zig()

  end subroutine inicializa

  subroutine finaliza (this)
    class (Monte_Carlo) :: this

     ! Libera memoria
    deallocate(this % R)
 end subroutine finaliza
 

  subroutine metropolis_oo (this)
    class (Monte_Carlo) :: this
    type(Lenard_Jones):: l_j

    real(dp) :: rx, ry, rz
    real(dp) :: n_energy, deltaE
    real(dp):: pr, beta
    integer :: iPart, i
    real(dp), dimension(3)  :: ri_vec

    beta = 1.0 / (this % params % gT * K_B_KJ)
    !beta = 1.0 / (gT * kb)

    allocate(this % abs_ener( this % params % gNtime))

    n_energy = l_j % potencial() 

        do i=1, this % params % gNtime
            iPart = rand_int(this % params % gNpart)
            ! guardo el valor original de la energia
            rx = this % R(1, iPart) 
            ry = this % R(2, iPart)
            rz = this % R(3, iPart)

            ! actualizo el arreglo de posiciones con una nueva posicion
            this % R(1, iPart) = this % R(1, iPart) + .2*this % params % gSigma * (uni() - 0.5)
            this % R(2, iPart) = this % R(2, iPart) + .2*this % params % gSigma * (uni() - 0.5)
            this % R(3, iPart) = this % R(3, iPart) + .2*this % params % gSigma * (uni() - 0.5)
     
            ! llamo a condiciones períodicas de contorno
            call cpc(iPart)

            ! calculo de la variacion de energia
            !deltaE = delta_poten_lj(iPart, rx, ry, rz)
            ri_vec = (/rx, ry, rz/)
            deltaE = l_j % delta_potencial(iPart, ri_vec)

            if (deltaE .lt. 0) then
                n_energy = n_energy + deltaE 
            else
                pr = exp(-beta * deltaE)
                if (uni() .lt. pr) then
                    n_energy = n_energy + deltaE    
                else
                    this % R(1, iPart) = rx
                    this % R(2, iPart) = ry
                    this % R(3, iPart) = rz
                endif        
            endif        
            this % abs_ener(i) = n_energy
            print *, n_energy 
        enddo

  endsubroutine metropolis_oo

endmodule mc
