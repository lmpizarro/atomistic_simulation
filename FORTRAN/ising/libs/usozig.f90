module usozig

  use ziggurat

  implicit none

contains

!==============================================================================
! Inicializa el generador de numeros aleatorios
!==============================================================================

  subroutine inic_zig()
    
    logical :: es
    integer :: seed   

    inquire(file='seed.dat',exist=es)
      if(es) then
        open(unit=10,file='seed.dat',status='old')
        read(10,*) seed
        close(10)
        print *,"  * Leyendo semilla de archivo seed.dat"
      else
        seed = 24583490
      end if

    call zigset(seed)

  end subroutine inic_zig

!==============================================================================
! Escribir la última semilla para continuar con la cadena de numeros aleatorios
!==============================================================================

  subroutine fin_zig()

    integer :: seed

    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
 
  end subroutine fin_zig

!==============================================================================
! Da un vector de numeros uniformemente distribuidos en [0,1] 
!==============================================================================

  function uni_vec(N) result(un)
  
    integer, intent(in) :: N
    real, dimension(N) :: un
    integer :: j

    do j=1,N
      un(j) = uni() 
    end do
  end function uni_vec

!==============================================================================
! Probabilidad uniforme tipo Bernoulli: enteros -1 o 1  
!==============================================================================

  function uni_2st() result(s)
    
    integer :: s
     
    if (uni()<= 0.5) then
      s = -1
    else
      s = 1
    end if

  end function uni_2st

!==============================================================================
! Probabilidad uniforme tipo Bernoulli: enteros -1 o 1  
!==============================================================================

  function rand_int(N) result(i_sel)
    ! Función que selecciona a un número al azar entre 1, ..., N posibilidades
    ! Asume que todos son equiprobables
      integer, intent(in)    :: N     ! Cantidad total 
      integer                :: i_sel ! Entero seleccionado

      i_sel = floor(uni()*N + 1)
  
  end function rand_int


end module usozig

