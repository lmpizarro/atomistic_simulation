module param_md 
  use types, only: dp
  implicit none

  real(dp),  dimension(:,:), allocatable  :: gR    ! Posicion de las partículas

  real(dp),  dimension(:,:), allocatable  :: gCombSigma    ! Matriz con sigma
  real(dp),  dimension(:,:), allocatable  :: gCombEpsilon  ! Matriz con epsilon
  real(dp), allocatable :: gRc2(:,:)  ! Cuadrado del radio de corte
  real(dp), allocatable :: gPot_Cut(:,:)  ! Potencial L-J en el radio de corte
  real(dp),  dimension(:), allocatable    :: gMasa         ! Vector de masas
  integer, allocatable :: gIndice_elemento(:)      ! Indice de partícula


  type Parametros
    real(dp) :: gTemperatura   ! Temperatura del sistema
    real(dp) :: gDt            ! Paso de tiempo de la corrida
    integer  :: gNtime         ! Cantidad de pasos de tiempo
    integer  :: gNmed          ! Cantidad de pasos entre mediciones
    real(dp) :: gGamma         ! Parámetro para el termostato de Langevin  
    integer  :: gNespecies      ! Cantidad de especies en el sistema
    integer  :: gNpart         ! Cantidad de partículas del sistema
    integer  :: gNperx         ! Cantidad de partículas del sistema
    integer  :: gNpery         ! Cantidad de partículas del sistema
    integer  :: gNperz         ! Cantidad de partículas del sistema
    real(dp) :: gParam_red
  contains
    procedure :: leer => read_parameters
  end type Parametros
contains

  subroutine read_parameters(this)
    class (Parametros) :: this
    logical :: es
    integer :: i, j
 
    inquire(file='./config.par',exist=es)
    if(es) then
      open(unit=10,file='./config.par',status='old')
      !lectura de los parámetros de corrida
      read(10,*) this%gTemperatura, this%gDt, this%gNtime, this%gNmed,&
         this%gGamma
      !lee cantidad de especies
      read(10,*) this%gNespecies
      
      ! alloca los arrays necesarios para leer configuraciones
      allocate(gCombSigma(1:this%gNespecies, 1:this%gNespecies))
      allocate(gCombEpsilon(1:this%gNespecies, 1:this%gNespecies))
      allocate(gRc2(1:this%gNespecies, 1:this%gNespecies))
      allocate(gPot_Cut(1:this%gNespecies, 1:this%gNespecies))
      allocate(gMasa(1:this%gNespecies))

      do i = 1, this%gNespecies
        read (10,*) (gCombSigma (i,j), j = 1,this%gNespecies)
      enddo

      do i = 1, this%gNespecies
        read (10,*) (gCombEpsilon (i,j), j = 1,this%gNespecies)
      enddo

      do i = 1, this%gNespecies
        read (10,*) (gRc2 (i,j), j = 1,this%gNespecies)
      enddo

      do i = 1, this%gNespecies
        read (10,*) (gPot_Cut (i,j), j = 1,this%gNespecies)
      enddo

      do i = 1, this%gNespecies
        read (10,*) gMasa (i)
      enddo

      read(10,*) this%gNpart, this%gNperx, this%gNpery, this%gNperz,&
        this%gParam_red

      ! alloca para leer posiciones de partículas
      allocate(gR(1:3, 1:this%gNpart))
      allocate(gIndice_elemento(1:this%gNpart))

      do i = 1, this%gNpart
        read (10,*) (gR (i,j), j = 1,3), gIndice_elemento(i)
      enddo

      close(10)
    else
       stop 1212      
    end if

    write(*, 100) this%gTemperatura, this%gDt, this%gNtime, this%gNmed,&
         this%gGamma
    write(*, 200) this%gNespecies

    do i = 1, this%gNespecies
      write (*, 700) (gCombSigma (i,j), j = 1,this%gNespecies)
    enddo

    do i = 1, this%gNespecies
      write (*, 700) (gCombEpsilon (i,j), j = 1,this%gNespecies)
    enddo

    do i = 1, this%gNespecies
      write (*, 700) (gRc2 (i,j), j = 1,this%gNespecies)
    enddo

    do i = 1, this%gnespecies
      write (*, 700) (gpot_cut (i,j), j = 1,this%gnespecies)
    enddo

    do i = 1, this%gnespecies
      write (*, 300) (gMasa (i))
    enddo

    write(*,200) this%gNpart, this%gNperx, this%gNpery, this%gNperz,&
        this%gParam_red

    do i = 1, this%gNpart
      write (*,800) (gR (i,j), j = 1,3), gIndice_elemento(i)
    enddo


    100 format (F14.8,F14.8,I10,I10,F14.8)
    200 format (I10, I10, I10, I10, F14.8)
    300 format (F14.8)
    700 format (100g15.5)
    800 format (3g15.5, I10)
  endsubroutine   

endmodule param_md
