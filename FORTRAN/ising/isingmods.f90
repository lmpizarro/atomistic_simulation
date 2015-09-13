module isingmods 

  use globales,            only: dp
  use ziggurat,            only: uni
  use usozig,              only: inic_zig, uni_2st, rand_int, fin_zig   
  use io_parametros,       only: escribe_estado, lee_estado

  implicit none

  private

!===============================================================================
! VARIABLES Y CONSTANTES  
!===============================================================================

  ! Variables
  integer, dimension(:,:), allocatable  :: RED      ! Matriz de spines
  real(dp)                              :: Eng      ! Energía del sistema
  real(dp)                              :: Mag      ! Magnetización
  integer                               :: N_R      ! Columnas de RED
  integer                               :: M_R      ! Filas de RED
  real(dp)                              :: Tem      ! Temperatura
  real(dp)                              :: beta     ! beta = 1/kT
  real(dp)                              :: Jac      ! Constante de acoplamiento
  integer                               :: K_tot    ! Cantidad de ciclos
  ! Constantes
  real(dp), parameter                   :: k_b= 1_dp ! Constante de Boltzman   
  real(dp), parameter                   :: mu = 1_dp ! Momento magnético

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
        Tem = 1.0_dp
        Jac = 1_dp
        K_tot = 100
     end if

    print *,"* Parámetros utilizados: ", N_R, M_R, Tem, Jac, K_tot

  end subroutine read_parameters

!===============================================================================
! INICIALIZA PARAMETROS
!===============================================================================

  subroutine inicializacion()
  
    integer  :: i,j
    logical  :: es

    ! Inicializa al generador de números aleatorios
    call inic_zig()
 
    ! Inicializa la matriz de spines 
    ! Define los indices pensando en las condiciones periodicas de contorno
    allocate( RED(0:N_R+1,0:M_R+1) )
    ! Trata de leer el estado inicial
    call lee_estado(es,RED,N_R,M_R)
    ! Si no leyó el estado, lo crea
    if (.not. es) then 
      do j = 1, M_R              ! Recordar que es column-major order
        do i = 1, N_R
         ! RED(i,j) = uni_2st()  ! Spines aleatorios
          RED(i,j) = 1           ! Todos los spines para arriba
        end do
      end do
    end if

    ! Aplica condiciones periódicas de contorno
    call cond_contorno(RED)  
    
    ! Inicializa beta
    beta = 1.0_dp/(k_b*Tem)

  end subroutine inicializacion 

!===============================================================================
! APLICA CONDICIONES PERIODICAS DE CONTORNO A LA MATRIZ 
!===============================================================================

  subroutine cond_contorno(A)

    integer, intent(inout), dimension(0:N_R+1,0:M_R+1) :: A  ! Matriz de spines 

    A(0,:)     = A(N_R,:)   ! Fila 0 igual a fila M_R
    A(N_R+1,:) = A(1,:)     ! Fila N_R+1 igual a fila 1 
    A(:,0)     = A(:,M_R)   ! Columna 0 igual a columna M_R
    A(:,M_R+1) = A(:,1)     ! Columna M_R+1 igual a columna 1

  end subroutine cond_contorno 

!===============================================================================
! CALCULA LA ENERGIA Y LA MAGNETIZACION
!===============================================================================

  subroutine calcula_EM() 

    real(dp)           :: E       ! Energía 
    real(dp)           :: M       ! Magnetización
    integer            :: i, j
    
    ! Inicializo las constantes. No conviene hacerlo al declararlas.
    E = 0.0_dp
    M = 0.0_dp

    do j = 1, M_R
      do i = 1, N_R
        ! Calcula la energía
        E = E - 0.5_dp*Jac*RED(i,j)*( RED(i+1,j) + RED(i-1,j) + RED(i,j+1) + RED(i,j-1) )
        ! Calcula la magnetización (sin el mu)
        M = M + RED(i,j)               
      end do
    end do 
    ! Escribo la energía total
    Eng = E
    ! Corrijo para obtener la magnetización
    Mag = mu*M             
  
    print *, '* Valores Iniciales: E = ', Eng, 'M = ', Mag
 
  end subroutine calcula_EM

!===============================================================================
! ALGORITMO DE METROPOLIS 
!===============================================================================

  subroutine metropolis()

    real(dp)     :: E_k    ! Energía del nuevo estado
    integer      :: k, i, j
    logical      :: acept  ! Flag para saber cuándo se acepta un estado

    
    ! Abro archivo para escribir los datos
    open(unit=20,file='energia.dat',status='unknown')
    open(unit=30,file='magneti.dat',status='unknown')

    do k = 1, K_tot

      ! Se genera el nuevo estado invirtiendo un spin al azar
      i = rand_int(N_R)     ! Genera al azar un entero para la fila
      j = rand_int(M_R)     ! Genera al azar un entero para la columna
     
      ! Calculo la energía del nuevov estado 
      E_k = Eng + 2.0_dp*Jac*RED(i,j)* ( RED(i-1,j) + RED(i+1,j) + RED(i,j-1) + RED(i,j+1) )
      
      ! ----------------- ACEPTACION-RECHAZO -----------------------------------
      ! Calculo la energía y la magnetización del nuevo estado
      ! Condiciones de aceptación-rechazo
      acept = .FALSE.                      

      if (E_k < Eng) then                  ! Si el nuevo estado es menos energético
        ! Acepto el nuevo estado, actualizo variables
        acept = .TRUE.
      else                                 ! Si el nuevo estado es más energético
        ! Acepto el estado con probabilidad e^{-\beta \Delta E}
        if ( uni() < exp(-beta*(E_k-Eng)) )  then
            acept = .TRUE.
        end if
      end if
      !------------------------------------------------------------------------

      ! Aplico la condición de contorno 
      if ( acept .eqv. .TRUE.) then      ! Si es aceptado
        RED(i,j) = -RED(i,j)             
        Eng = E_k
        Mag = Mag + 2.0*mu*RED(i,j) 
        if ( i==1 .or. i==N_R .or. j==1 .or. j==M_R) then ! Si es spin periférico
          call cond_contorno(RED)
        end if
      end if

      ! Guardo los datos de k, Eng y Mag
      write(20,100)  Eng
      write(30,100)  Mag
  100 format(F8.1)

    end do 

    ! Cierro archivo de salida
    close(20)   
    close(30)   

    ! Trata de escribir el último estado a un archivo
    call escribe_estado(RED)

    ! Informa el estado final
    print *, '* Valores Finales: E = ', Eng, 'M = ', Mag

  end subroutine metropolis 

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
