module estadistica
  
  use types,               only: dp

  implicit none
 
  private

  public   :: vector_mean_var_dir, histograma_vec

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
