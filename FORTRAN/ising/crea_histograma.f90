program main

!===============================================================================
! PROGRAMA PARA CONSTRUIR HISTOGRAMAS
!===============================================================================
!
! Se ingresan los parámetros a través de la línea de comandos
! El uso es:
!    -i <nombre del archivo de entrada> (obligatorio)
!    -n <cantidad de bines>   (opcional, por defecto es 100)
!    -l <limite inferior de los bines> (opcional)
!    -u <limite superior de los bines> (opcional)
!         Si no se incluyen los -l y -u se los toma como el
!         mínimo y el máximo de los datos ingresados

use globales,            only: dp
use estadistica,         only: histograma
use io_parametros,       only: read_command_line
use strings,             only: str_to_int, str_to_real

implicit none

integer                              :: N
real(dp),dimension(:), allocatable   :: cuentas
real(dp),dimension(:), allocatable   :: bins
real(dp),dimension(2)                :: lim_op
integer                              :: i

character(100) :: nombre, N_lin_com, bin_min, bin_max

lim_op = [1_dp,5_dp]

! Lee el archivo de entrada de la línea de comando
call read_command_line(nombre,'-i')

if (nombre == '') then
  write(*,*) 'Error al ingresar el nombre del archivo'
  stop
end if

! Lee la cantidad de bines
call read_command_line(N_lin_com,'-n')

if (N_lin_com /='') then
  N = str_to_int(trim(N_lin_com))
else
  N = 100
end if

call read_command_line(bin_min,'-d')
call read_command_line(bin_max,'-u')

if (bin_min /='') then
  lim_op(1) = str_to_real(trim(bin_min))
end if

if (bin_max /='') then
  lim_op(2) = str_to_real(trim(bin_max))
end if

allocate(cuentas(N),bins(N))

if ( (bin_max =='') .or. (bin_max == '') ) then
  call histograma(cuentas,bins,nombre,N)
else
  call histograma(cuentas,bins,nombre,N,lim_op)
end if

write(*,'(2(E12.5))') ( bins(i), cuentas(i), i=1,N)

deallocate(cuentas,bins)

end program main 

