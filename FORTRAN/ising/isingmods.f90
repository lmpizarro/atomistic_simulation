module isingmods 

  use usozig
  use globals

  implicit none

  private
!===============================================================================
! VARIABLES Y CONSTANTES  
!===============================================================================
  ! Variables
  integer                              :: N_R      ! Columnas de RED
  integer                              :: M_R      ! Filas de RED
  real                                 :: Tem      ! Temperatura
  real                                 :: beta     ! beta = 1/kT
  real                                 :: Jac      ! Constante de acoplamiento
  integer                              :: K_tot    ! Cantidad de ciclos
  ! Constantes
  real, parameter           :: k_b = 1.            ! Constante de Boltzman   
  real, parameter           :: mu = 1.             ! Momento magnético

  public :: read_parameters, inicializacion, calcula_EM, metropolis, finalizacion

contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
  
    logical :: es
     
    inquire(file='parametros.dat',exist=es)
      if(es) then
        open(unit=10,file='parametros.dat',status='old')
        read(10,*) N_R, M_R, Tem, Jac, K_tot
        close(10)
      else
        N_R = 10
        M_R = 10
        Tem = 1.0
        Jac = 1
        K_tot = 100
     end if

    print *,"* Parámetros utilizados: ", N_R, M_R, Tem, Jac, K_tot

  end subroutine read_parameters

!===============================================================================
! INICIALIZA PARAMETROS
!===============================================================================

  subroutine inicializacion()
  
    integer  :: i,j

    ! Inicializa al generador de números aleatorios
    call inic_zig()
 
    ! Inicializa la matriz de spines 
    ! Define los indices pensando en las condiciones periodicas de contorno
    allocate( RED(0:N_R+1,0:M_R+1) )

    do j = 1, M_R               ! Recordar que es column-major order
      do i = 1, N_R
        RED(i,j) = uni_2st()    ! Spines aleatorios
       ! RED(i,j) = 1            ! Todos los spines para arriba
      end do
    end do

    ! Aplica condiciones periódicas de contorno
    call cond_contorno(RED)  

    ! Inicializa beta
    beta = 1/(k_b*Tem)

  end subroutine inicializacion 

!===============================================================================
! APLICA CONDICIONES PERIODICAS DE CONTORNO A LA MATRIZ 
!===============================================================================

  subroutine cond_contorno(A)

    integer, intent(inout), dimension(0:N_R+1,0:M_R+1) :: A  ! Matriz de spines 

    A(0,0:M_R+1)     = A(N_R,0:M_R+1)   ! Columna 0 igual a columna M_R
    A(N_R+1,0:M_R+1) = A(1,0:M_R+1)     ! Columna N_R+1 igual a columna 1 
    A(0:N_R+1,0)     = A(0:N_R+1,M_R)   ! Fila 0 igual a fila M_R
    A(0:N_R+1,M_R+1) = A(0:N_R+1,1)     ! Fila M_R+1 igual a fila 1

  end subroutine cond_contorno 

!===============================================================================
! CALCULA LA ENERGIA Y LA MAGNETIZACION
!===============================================================================

  subroutine calcula_EM(E,M,A)

    real, intent(out)                               :: E    ! Energía 
    real, intent(out)                               :: M    ! Magnetización
    integer, intent(in), dimension(0:N_R+1,0:M_R+1) :: A    ! Matriz de spines  
  
    integer   :: i, j
    
    ! Inicializo las constantes. No conviene hacerlo al declararlas.
    E = 0.
    M = 0.

    do j = 1, M_R
      do i = 1, N_R
        ! Calcula la energía
        E = E - 0.5*Jac*A(i,j)*( A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1) )
        ! Calcula la magnetización (sin el mu)
        M = M + A(i,j)               
      end do
    end do 
    ! Corrijo para obtener la magnetización
    M = mu*M            ! Revisar por las dudas 

  end subroutine calcula_EM

!===============================================================================
! ALGORITMO DE METROPOLIS 
!===============================================================================

  subroutine metropolis()
    ! Magnitudes temporales par los estados temporarios en el paso k-ésimo
    real                                 :: E_k    ! Energía 
    real                                 :: M_k    ! Magnetización
    integer,  dimension(0:N_R+1,0:M_R+1) :: RED_k  ! Matriz de spines  

    integer   :: k, i, j

    ! Inicializo las constantes. No conviene hacerlo al declararlas.
    E_k = 0.
    M_k = 0.

    do k = 1, K_tot
      ! Copia el nuevo estado temporal
      RED_k = RED     
      ! Se genera el nuevo estado invirtiendo un spin al azar
      i = rand_int(N_R)     ! Genera al azar un entero para la columna
      j = rand_int(M_R)     ! Genera al azar un entero para la fila
      ! Se da vuelta al spin seleccionado
      RED_k(i,j) = -RED(i,j)
      ! Aplico la condición de contorno sólo si cambio un spin periférico
      if ( i==1 .or. i==N_R .or. j==1 .or. j==M_R) then
        call cond_contorno(RED_k)
      end if  
      ! ----------------- ACEPTACION-RECHAZO -----------------------------------
      ! Calculo la energía y la magnetización del nuevo estado
      call calcula_EM(E_k,M_k,RED_k)
      ! Condiciones de aceptación-rechazo
      if (E_k < Eng) then       ! Si el nuevo estado es menos energético
        ! Acepto el nuevo estado, actualizo variables
        RED = RED_k
        Eng = E_k
        Mag = M_k 
      else                      ! Si el nuevo estado es más energético
        ! Acepto el estado con probabilidad e^{-\beta \Delta E}
        if ( uni() < exp(-beta*(E_k-Eng)) )  then
          RED = RED_k           ! Acepto el estado
          Eng = E_k
          Mag = M_k
        end if
      end if
      ! Guardo los datos de k, Eng y Mag
      print *, k, Eng, Mag
    end do 


  end subroutine metropolis 

!================================================================================
! ESCRIBE LOS RESULTADOS A UN ARVHIVO
!================================================================================ 
  
  subroutine escribe_datos(nombre,x,y)

    character (len=*), intent(in) :: nombre
    real, intent(in), dimension(:) :: x
    real, intent(in), optional, dimension(:) :: y
    integer :: j, N
    
    N = size(x)
    if (present(y))  then
      open(unit=40,file=nombre,status='unknown')
      write(40,'(2F15.5)') ([x(j), y(j)] ,j=1,N)
      close(40)
    else
      open(unit=40,file='x_'//nombre,status='unknown')
      write(40,'(F15.5)') (x(j) ,j=1,N)
      close(40) 
    end if

  end subroutine escribe_datos

!===============================================================================
! FINALIZA PARAMETROS
!===============================================================================

  subroutine finalizacion()
  
    ! Libera memoria
    deallocate(RED)

    ! Escribe la semmilla del generador de números aleatorios
    call fin_zig()    

  end subroutine finalizacion

end module isingmods 
