module io_parametros

  use globales,     only: dp
  use strings,      only: starts_with

  implicit none

contains

!================================================================================
! LEE EL ESTADO DE SPINES EN UN ARCHIVO 
!================================================================================ 
  
  subroutine lee_estado(es,A,N,M)

  integer, intent(out), dimension(:,:), allocatable ::  A  ! Matriz de spines
  logical, intent(out)  :: es
  integer, intent(in)   :: N,M

  allocate(A(0:N+1,0:M+1)) 
   
    inquire(file='ultimo_estado.dat',exist=es)
    if (es) then
      open(unit=10,file='ultimo_estado.dat',status='old')
      read(10,*) A
      close(10)
      print *, '* Se lee estado inicial desde el archivo <ultimo_estado.dat>'
    end if
  
  end subroutine lee_estado


!================================================================================
! ESCRIBE EL ESTADO DE SPINES EN UN ARCHIVO 
!================================================================================ 
  
  subroutine escribe_estado(A)

    logical    :: es
    integer, intent(in), dimension(:,:) ::  A  ! Matriz de spines
     
    inquire(file='ultimo_estado.dat',exist=es)
    if ( .not. es) then
      open(unit=10,file='ultimo_estado.dat',status='new')
      write(10,*) A
      close(10)
      print *, '* Se escribe el último estado en el archivo <ultimo_estado.dat>'
    end if

  end subroutine escribe_estado

!================================================================================
! ESCRIBE EL ESTADO DE SPINES EN UN ARCHIVO 
!================================================================================ 
  
  subroutine escribe_aceptaciones(T,N,N_tot)

    real(dp), intent(in)      :: T         ! Temperatura
    integer, intent(in)       :: N, N_tot  ! Número de aceptaciones y totales
     
      open(unit=10,file='aceptaciones.dat',status='unknown')
      write(10,'(F6.2,4X,I10,4X,I10)') T, N, N_tot
      close(10)
      print *, '* Se escribe el número de aceptaciones en el archivo <aceptaciones.dat>'

  end subroutine escribe_aceptaciones

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line(salida,tipo)

    character(*), intent(out) :: salida   ! Valor del argumento buscado
    character(*), intent(in)  :: tipo     ! Identificación del argumento buscado
    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    character(150), allocatable :: argv(:) ! command line arguments

    ! Check number of command line arguments and allocate argv
    argc = COMMAND_ARGUMENT_COUNT()

    ! Allocate and retrieve command arguments
    allocate(argv(argc))
    do i = 1, argc
      call GET_COMMAND_ARGUMENT(i, argv(i))
    end do

    ! Process command arguments
    salida = ''
    i = 1
    do while (i <= argc)
      ! Check for flags
      if (starts_with(argv(i), "-")) then
        if (tipo==argv(i)) then
          i = i + 1
          salida = argv(i)
        end if
      end if
      ! Increment counter
      i = i + 1
    end do

    ! Free memory from argv
    deallocate(argv)

  end subroutine read_command_line

end module io_parametros
