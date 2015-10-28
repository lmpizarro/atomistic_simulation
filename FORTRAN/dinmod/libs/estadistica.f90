module estadistica
  
  use globales,            only: dp

  implicit none

contains

!===============================================================================
! CALCULA LOS DATOS PARA GRAFICAR UN HISTOGRAMA - DE FORMA VECTORIAL
!===============================================================================

  subroutine histograma_vec(cuentas,bins,x,N,lim_op)

    real(dp), dimension(1:), intent(out)     :: cuentas! Cuentas normalizadas
    real(dp), dimension(1:), intent(out)     :: bins   ! bines utilizados
    real(dp), dimension(1:), intent(in)      :: x      ! Datos a binnear
    integer, intent(in)                      :: N      ! Cantidad de datos
    real(dp), dimension(2),  optional        :: lim_op ! [lim_inf, lim_sup]
    integer, dimension(:), allocatable       :: Ncuen  ! Cuentas sin normalizar

    real(dp), dimension(2)                   :: lim    ! Límites de los bines
    real(dp)                                 :: dbin   ! Ancho del bin
    real(dp), dimension(N+1)                 :: rango  ! Comienzo de cada bin
    integer                                  :: j
 
    if (present(lim_op) .neqv. .TRUE.) then 
      lim = [minval(x), maxval(x)] 
    else
      lim = lim_op
    end if
    
    dbin = (lim(2)-lim(1))/N
    ! Originalmente estaba float(j)
    rango = [ (lim(1)+real(j,dp)*dbin,j=0,N) ] 

    allocate(Ncuen(N))
    
    Ncuen=0
    ! Todas los puntos menores al  límite inferior van al primer bin
    Ncuen(1) = count( ( x<rango(2) ) )   
    do j=2,N
      Ncuen(j) = count( ( x<rango(j+1) ) .and. ( x>=rango(j) ) )
    end do
   ! Todslos los puntos mayores al límite superior van al último bin
    Ncuen(N) = Ncuen(N) + count( x>=rango(N+1) )   
   ! Normalizo las cuentas
    cuentas = real(Ncuen(1:N),dp) / sum(Ncuen) / dbin 
   ! Uso un punto centrado en el bin para graficar 
    bins = rango(1:N) + dbin/2    
    return

    deallocate(Ncuen)  

  end subroutine histograma_vec

!===============================================================================
! CALCULA LOS DATOS PARA GRAFICAR UN HISTOGRAMA - PUNTO POR PUNTO
!===============================================================================

  subroutine histograma(cuentas,bins,nombre,N,lim_op)

    real(dp), dimension(1:N), intent(out)    :: cuentas! Cuentas normalizadas
    real(dp), dimension(1:N), intent(out)    :: bins   ! bines utilizados
    character (len=*), intent(in)            :: nombre ! Nombre del archvo a leer
    integer, intent(in)                      :: N      ! Cantidad de datos
    real(dp), dimension(2),  optional        :: lim_op ! [lim_inf, lim_sup]

    integer, dimension(1:N)                  :: Ncuen  ! Cuentas sin normalizar
    real(dp), dimension(2)                   :: lim    ! Límites de los bines
    real(dp)                                 :: dbin   ! Ancho del bin
    real(dp), dimension(0:N)                 :: rango  ! Comienzo de cada bin
    real(dp)                                 :: x      ! Dato leidos
    integer                                  :: estado ! Estado del valor leido
    integer                                  :: j,k

    if (present(lim_op) .neqv. .TRUE.) then 
      ! Si no se dan los límiter, se lee el archivo para detectar los límites
      ! mínimos y máximos de los datos
      call primer_leida(lim,nombre)  
    else
      lim = lim_op
    end if

    ! Se define el ancho del bin
    dbin = (lim(2)-lim(1))/N
    ! Vector con los extremos de cada bin 
    ! Si se tienen N bines, 'rango' tendrá N+1 elementos
    rango = [ (lim(1)+real(j,dp)*dbin,j=0,N) ] 

    ! Abre el archivo
    open(20, file = nombre)
    ! Inicializa el vector donde se acumularán las cuentas
    Ncuen = 0
    do    ! Loop sobre todos los puntos que tiene el archivo
      read(20,*,iostat=estado) x
      if (estado > 0) then
        print *, 'ERROR al leer los datos'
        exit
      else if (estado < 0) then
        !print*, 'Lectura del archivo completa'
        exit
      else
        do k = 1, N     ! Loop sobre los valores de rango
          if ( x<rango(k) ) then
            ! Cuenta en el bin correspondiente y sale del loop
            Ncuen(k) = Ncuen(k) + 1
            exit
          end if
        end do
        ! Los mayores a lim(2) van en el ultimo bin
        if ( x>= rango(N) ) Ncuen(N) = Ncuen(N) + 1
      end if
    end do

    close(20)

    ! Normalizo las cuentas
    cuentas = real(Ncuen(1:N),dp) / sum(Ncuen) / dbin    
    ! Uso un punto centrado en el bin para graficar 
    bins = rango(1:N) - dbin/2    
    return

  contains

    subroutine primer_leida(limites,nombre)
    ! En el caso de no especificar los línmites, esta subrutina lee el archivo
    ! y busca los valores mínimos y máximo de los datos
      real(dp), dimension(2), intent(out)     :: limites ! [lim_inf, lim_sup]
      character (len=*), intent(in)           :: nombre  ! Nombre del archvo a leer

      integer                                 :: estado  ! Estado del valor leido
      real(dp)                                :: x, x_min , x_max 

      x_min = huge(real(dp))
      x_max = -huge(real(dp))

      open(20, file = nombre)

      do
        read(20,*,iostat=estado) x
        if (estado > 0) then
          print *, 'ERROR al leer los datos durante la determinación de los límites'
          exit
        else if (estado < 0) then
          limites = [x_min, x_max]
          !print*, 'Se determinaron correctamente los límites del histograma por defecto'
          exit
        else
          if (x < x_min) then
            x_min = x
          else if (x > x_max) then
            x_max = x
          end if
        end if
      end do

      close(20)

    end subroutine primer_leida 

  end subroutine histograma


!===============================================================================
! PROMEDIO Y DESVIACION ESTANDAR DE FORMA ITERATIVA
!===============================================================================
! Se evita cargar en memoria todos los datos a la vez. Se leen de a uno y se 
! calculan el promedio y la varianza normalizada con sqrt(1-N) recursivamente
! ADVERTENCIA: El método usado tiene problemas de redondeo para N muy grande.
! Se tendria que modificar el algoritmo pararealizar el promedio (i.e de a pares)
! Sin embargo, trabajando con doble precision se mejora bastante sin hacer nada.

  subroutine running_mean_var(x_mean,x_var,nombre)

    character (len=*), intent(in)  :: nombre          ! Nombre del archvo a leer
    real(dp), intent(out)          :: x_mean, x_var   ! Valor medio y varianza
    real(dp)                       :: x               ! Valor leido
    real(dp)                       :: x_mean_old      ! Valor temporal
    integer                        :: estado          ! Estado del valor leido
    integer                        :: N               ! Cantidad de datos leidos

    ! Abre el archivo
    open(20, file = nombre)
    ! Lee el primer punto e inicializa los valores
    read(20,*,iostat = estado) x
    x_mean = x
    x_var = -1.0_dp
    N =  1
    ! Loop que lee linea por linea al archivo hasta que encuentre el EOF
    do
      read(20,*,iostat=estado) x
      if (estado > 0) then
        print *, 'ERROR al leer los datos'
        exit
      else if (estado < 0) then
        print*, 'Lectura del archivo completa'
        exit
      else
      ! Guarda el valor medio viejo para usarlo cundo calcula la varianza
      x_mean_old = x_mean
      x_mean =  (N*x_mean + x)  / ( N + 1 )
      x_var = ( 1.0_dp - real(1,dp)/N )*x_var + ( N + 1 )*( x_mean - x_mean_old )**2
      N = N + 1
      end if
    end do

    close(20)

  end subroutine running_mean_var

!===============================================================================
! PROMEDIO Y DESVIACION ESTANDAR VECTORIZADO
!===============================================================================
! La desventaja es que debe cargar al archivo completo en memoria. La ventaja es
! que tiene menos errores de redondeo.

  subroutine vector_mean_var(x_mean,x_var,nombre)

    character (len=*), intent(in)       :: nombre        ! Nombre del archvo a leer
    real(dp), intent(out)               :: x_mean, x_var ! Valor medio y varianza
    real(dp)                            :: x             ! Valor leido
    real(dp), dimension(:), allocatable :: x_vec         ! Array para guardar datos
    integer                             :: N,i,estado 
  
    ! Esta parte solo esta para saber cuantos datos hay en el archivo
    ! Se podria evitar guardando el valor de N al comienzo.
    N = 0
    open(20, file = nombre)
    do
      read(20,*,iostat=estado) x
      if (estado > 0) then
        print *, 'ERROR al leer los datos'
        exit
      else if (estado < 0) then
        print*, 'Lectura del archivo completa'
        exit
      else
        N = N + 1         ! Cuenta los datos leidos
      end if
    end do
    close(20)
    ! Una vez conocido N se puede inicializar al vector donde se leeran los datos
    allocate(x_vec(N))
    ! Se vuelve a abrir el archivo y a leer los datos 
    open(20,file =nombre)
    read(20,*) (x_vec(i),i=1,N)
    close(20)
    ! Realiza las operaciones estadisticas de forma vectorial 
    x_mean=sum(x_vec)/N
    x_var= sum((x_vec - x_mean)**2) / (N-1)

    deallocate(x_vec)

  end subroutine vector_mean_var


!===============================================================================
! PROMEDIO Y DESVIACION ESTANDAR VECTORIZADO DIRECTO
!===============================================================================
! La desventaja es que debe cargar al archivo completo en memoria. La ventaja es
! que tiene menos errores de redondeo.
! Es igual a vector_mean_var() solo que opera sobre el vector de entrada
! y no sobre el archivo con los datos
  subroutine vector_mean_var_dir(x_mean,x_var,x)

    real(dp), intent(in), dimension(:)  :: x             ! Vector con datos entrada
    real(dp), intent(out)               :: x_mean, x_var ! Valor medio y varianza
    integer                             :: N 

    ! Tamaño del vector de entrada
    N = size(x,1)  

    ! Realiza las operaciones estadisticas de forma vectorial 
    x_mean=sum(x)/N
    x_var= sum((x - x_mean)**2) / (N-1)

  end subroutine vector_mean_var_dir

end module estadistica
