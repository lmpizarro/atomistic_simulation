module datos_problema

  use types, only: dp

  implicit none
  type Parametros
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
    real(dp) :: gSigma, gEpsil, gM

  contains
    procedure :: leer => read_parameters
  end type Parametros

contains

  subroutine read_parameters(this)
    class (Parametros) :: this
    logical :: es
 
    inquire(file='./parametros.dat',exist=es)
    if(es) then
      open(unit=10,file='./parametros.dat',status='old')
      read(10,*) this%gT, this%gNpart, this%gL, this%gDt, this%gNtime, this%gSigma,&
        this%gEpsil, this%gM, this%gNmed
      close(10)
    else
      print *, "Usando parametros por default"
      this%gT = 293.0_dp
      this%gNpart = 1000
      this%gL = 10.0_dp
      this%gDt = 0.1_dp
      this%gNtime = 100
      this%gSigma = 1.0_dp
      this%gEpsil = 1.0_dp
      this%gM = 1.0_dp
    end if
   
    write(*,'(a)') ''
    write(*,'(a)')      '************ PARAMETROS LEIDOS **************'
    write(*,'(a,F8.3)') '************ Temperatura           = ' , this%gT 
    write(*,'(a,I8)')   '************ Número de partículas  = ' , this%gNpart 
    write(*,'(a,F8.3)') '************ Masa de la partícula  = ' , this%gM
    write(*,'(a,F8.3)') '************ Lado del cubo         = ' , this%gL 
    write(*,'(a,F8.5)') '************ Paso temporal dt      = ' , this%gDt 
    write(*,'(a,I8)')   '************ Número de pasos dt    = ' , this%gNtime 
    write(*,'(a,I8)')   '************ Pasos mediciones      = ' , this%gNmed 
    write(*,'(a)')      '************ Potencial de L-J' 
    write(*,'(a,F8.4)') '************ Epsilon               = ' , this%gEpsil
    write(*,'(a,F8.4)') '************ Sigma                 = ' , this%gSigma
    write(*,'(a,F8.4)') '************ Masa                  = ' , this%gM
    write(*,'(a)')      '*********************************************'

 
  end subroutine read_parameters
end module datos_problema




