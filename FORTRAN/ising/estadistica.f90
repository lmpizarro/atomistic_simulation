module estadistica
  
  use globales,            only: dp

  implicit none

contains

!===============================================================================
! CALCULA LOS DATOS PARA GRAFICAR UN HISTOGRAMA 
!===============================================================================

  subroutine histograma(cuentas,bins,x,N,lim_op)

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

end module estadistica
