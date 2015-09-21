module isingmods 

  use globales,            only: dp
  use ziggurat,            only: uni
  use usozig,              only: inic_zig, uni_2st, rand_int, fin_zig   
  use io_parametros,       only: escribe_estado, lee_estado, read_command_line, &
                                 escribe_aceptaciones
  use strings,             only: str_to_int, str_to_real 
  use estadistica,         only: vector_mean_var

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

  public :: read_parameters, inicializacion, calcula_EM, metropolis, &
            finalizacion, hace_estadistica

contains

!===============================================================================
! LEE LOS DATOS DE ENTRADA DEl PROBLEMA 
!===============================================================================

  subroutine read_parameters()
  
    logical :: es
    character(100) :: K_tot_lin_com, Tem_lin_com
     
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
    
    ! Se fija si se ingresaron parametros por la linea de commandos
    ! El numero total de ciclos
    call read_command_line(K_tot_lin_com,'-k')
    if (K_tot_lin_com /='') K_tot = str_to_int(trim(K_tot_lin_com))
    ! La temperatura
    call read_command_line(Tem_lin_com,'-t')
    if (Tem_lin_com /='')  Tem   = str_to_real(trim(Tem_lin_com))
    print *, '------------- PARAMETROS UTILIZADOS -------------' 
    write(*, '(A,23X,I2,A,I2)') ' ---- Red de spines:' , N_R, ' x ', M_R
    write(*,'(A,25X,F6.3)')  ' ---- Temperatura: ', Tem
    write( *,'(A,10X,F5.3)') ' ---- Constante de acomplamiento J:', Jac
    write( *,'(A,1X,I10)')  ' ---- Número de pasos para Metrópolis::', K_tot
    print *, '-------------------------------------------------' 

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
    RED = 1
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
  
    write(*,'(A,F9.2,4X,A,F9.2)') ' * Valores Iniciales:  E =', Eng, 'M =', Mag

  end subroutine calcula_EM

!===============================================================================
! ALGORITMO DE METROPOLIS 
!===============================================================================

  subroutine metropolis()

    real(dp)     :: E_k    ! Energía del nuevo estado
    integer      :: k, i, j
    logical      :: acept  ! Flag para saber cuándo se acepta un estado
    integer      :: num_acept  ! Cantidad de puntos aceptados
    
    ! Inicializo el contador
    num_acept = 0
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

      ! En caso de ser aceptado, actualizo valores 
      if ( acept .eqv. .TRUE.) then    
        RED(i,j) = -RED(i,j)             
        Eng = E_k
        Mag = Mag + 2.0*mu*RED(i,j) 
        ! Aplico la condición de contorno sólo si cambió un spin periférico
        if ( i==1 .or. i==N_R .or. j==1 .or. j==M_R) then
          call cond_contorno(RED)
        end if
        num_acept = num_acept + 1
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
    ! Escribe la informacion sobre la cantdiad de estados aceptados
    call escribe_aceptaciones(Tem,num_acept,K_tot)

    ! Informa el estado final
    write(*,'(A,F9.2,4X,A,F9.2)') ' * Valores Finales:    E =', Eng, 'M =', Mag

  end subroutine metropolis 

!===============================================================================
! CALCULA VALOR MEDIO Y VARIANZA DE LOS ARCHIVOS DE SALIDA 
!===============================================================================

  subroutine hace_estadistica()

    real(dp)     :: E_media, M_media  ! Valor medio de la energía y la magnetización
    real(dp)     :: E_var , M_var     ! Varianza de la energía y la magnetización
  
    ! Calcula los valores medios y varianzas de la corrida
    ! Para la energía
    print *, '* Se abre archivo <energia.dat> para hacer estadistica'
    call vector_mean_var(E_media,E_var,'energia.dat')
  
    ! Para la magnetización
    print *, '* Se abre archivo <magneti.dat> para hacer estadistica'
    call vector_mean_var(M_media,M_var,'magneti.dat')

    ! Guarda los resultados en un archivo
    print *, '* Se guardan valores medios y varianza en <val_medios.dat>'
    open(unit=40,file='val_medios.dat',status='unknown')
    write(40,120) Tem, E_media, (E_var)*beta/Tem , M_media, (M_var)*beta
    120 format(F5.2, 4E15.7)
    close(40)

  end subroutine hace_estadistica

!===============================================================================
! FINALIZA PARAMETROS
!===============================================================================

  subroutine finalizacion()
      
    integer  :: i, j     ! contadores

    ! Escribe el ultimo estado en donde quedo el sistema
    ! En forma matricial (para verlo y graficarlo facilmente)
    open(unit=10,file='estado_mat.dat',status='unknown')
    do i= 0, N_R+1
       write(10,'(900I3)') (RED(i,j), j = 0, M_R+1)
    end do
    close(10)
   
    ! Escribe el ultimo estado  para que sea leido como estado inicial
    open(unit=20,file='estado.dat',status='unknown')
       write(20,*) RED
    close(20)

    ! Libera memoria
    deallocate(RED)

    ! Escribe la semmilla del generador de números aleatorios
    call fin_zig()    

  end subroutine finalizacion

end module isingmods 
