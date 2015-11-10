module procesamiento

  use types,             only: dp
  use globales,          only: gNpart, gL, gRho, gDt, gNtime, gNmed, gNpart, gT 
  use estadistica,       only: vector_mean_var_dir

  implicit none

  private

  public   :: hace_estadistica 

contains

  !===============================================================================
  ! ESTADISTICA SOBRE VECTORES Y ESCRITURA DE RESULTADOS  
  !===============================================================================
  ! Se calculan valores medios y desviacions estándares
  ! Se escriben los resultados al archivo <val_medios.dat> 

  subroutine hace_estadistica(x,y,z)
    
    real(dp), dimension(:), intent(in)      :: x, y  ! Datos en 1D para hacer estadística
    real(dp), dimension(:,:), intent(in)    :: z     ! Datos en 2D para hacer estadística

    real(dp)                         :: x_mean, y_mean ! Valor medio de los datos
    real(dp)                         :: x_std, y_std   ! Desviación estandar de los datos
    real(dp), dimension(1:size(z,1)) :: z_mean, z_std  ! para 2D
    integer                          :: i
    
    ! Calcula valores medios para el primer parámetro (presion) 
    write(*,*) '* Se hace estadistica con los valores de presión'
    call vector_mean_var_dir(x_mean,x_std,x)
    x_std = sqrt(x_std)

    ! Calcula valores medios para el segundo parámetro (temperatura) 
    write(*,*) '* Se hace estadistica con los valores de temperatura'
    call vector_mean_var_dir(y_mean,y_std,y)
    y_std = sqrt(y_std)
    
    ! Parche para calcular estadistica por columnas en una matriz
    ! Hacerlo más prolijo en algún momento
    write(*,*) '* Se hace estadistica con los valores de energías'
    z_mean = sum(z,2) / size(z,2) 

    do i = 1, size(z,1)
      z_std(i) =  sum( (z(i,:)-z_mean(i))**2 ) / (size(z,2) -1) 
    end do
    z_std = sqrt(z_std)
    ! Fin parchado 

    !-------------------------------------------------------------------------------------
    ! ESCRITURA DE DATOS
    !------------------------------------------------------------------------------------

    ! Se guardan los resultados en un archivo
    write(*,*) '* Se guardan valores medios y desv. estandar en <val_medios.dat>'

    open(unit=10, file='val_medios.dat', status='unknown')
    ! Primera parte con datos de la corrida   
    write(10,100)   'Temp', '#_part', 'Lado_cubo', 'Dens', 'dt', '#_pasos',  '#_pasos_med.'
    100 format (2X,7(A,8X))
    write(10,200) gT, gNpart, gL, gRho, gDt, float(gNtime), float(gNmed)
    200 format (F6.2,10X,I6,9X,F11.5,5X,F11.5,5X,3(E9.2,5X))
    write(10,*)

    ! Segunda linea con los valores medios y desviaciones estándars
    write(10,300) 'mean(p)', 'std(p)' , 'mean(T)', 'std(T)'
    300 format (4X,4(A,12X))
    write(10,400) x_mean, x_std, y_mean, y_std
    400 format (4(E14.7,4X))
    write(10,*)
 
    ! Tercera linea con los valores medios y desviaciones estandars
    write(10,500) 'mean(E_pot)', 'std(E_pot)', 'mean(E_kin)', 'std(E_kin)', 'mean(E_tot)', 'std(E_kin)'
    500 format (2X,6(A,8X))
    write(10,600) (  (/z_mean(i), z_std(i)/)  , i=1,3) 
    600 format (6(E14.7,4X))

    close(10)

  end subroutine hace_estadistica 

end module procesamiento
