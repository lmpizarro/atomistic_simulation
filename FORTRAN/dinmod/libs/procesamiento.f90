module procesamiento

  use types,             only: dp
  use globales,          only: gNpart, gL, gRho, gDt, gNtime, gNmed 
  use estadistica,       only: vector_mean_var_dir

  implicit none

  private

  public   :: hace_estadistica 

contains

  subroutine hace_estadistica(x,y)
    
    real(dp), dimension(:), intent(in)      :: x, y   ! Datos en 1D para hacer estadística

    real(dp)                :: x_mean, y_mean ! Valor medio de los datos
    real(dp)                :: x_std, y_std  ! Desviación estandar de los datos
    
    ! Calcula valores medios para el primer parámetro (presion) 
    write(*,*) '* Se hace estadistica con los valores de presión'
    call vector_mean_var_dir(x_mean,x_std,x)
    x_std = sqrt(x_std)

    ! Calcula valores medios para el segundo parámetro (temperatura) 
    write(*,*) '* Se hace estadistica con los valores de temperatura'
    call vector_mean_var_dir(y_mean,y_std,y)
    y_std = sqrt(y_std)

    ! Se guardan los resultados en un archivo
    write(*,*) '* Se guardan valores medios y desv. estandar en <val_medios.dat>'

    open(unit=10, file='val_medios.dat', status='unknown')
    ! Primera parte con datos de la corrida   
    write(10,100)   '# partículas', 'Lado cubo', 'Densidad', 'dt', '# pasos',  '# pasos med.'
    100 format (6(A,8X))
    write(10,200) gNpart, gL, gRho, gDt, float(gNtime), float(gNmed)
    200 format (2X,I6,10X,F11.5,4X,F11.5,4X,3(E9.2,6X))
    write(10,*)
    ! Segunda linea con los valores medios y desviaciones estándars
    write(10,300) 'mean(p)', 'std(p)' , 'mean(T)', 'std(T)'
    300 format (4X,4(A,12X))
    write(10,400) x_mean, x_std, y_mean, y_std
    400 format (4(E14.7,4X))

    close(10)

  end subroutine hace_estadistica 

end module procesamiento
