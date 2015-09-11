
module io_parametros
 
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

end module io_parametros
