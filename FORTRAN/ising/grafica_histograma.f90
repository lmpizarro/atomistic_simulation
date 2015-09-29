program main
!==========================================================
! PROGRAMA DE PRUEBA PARA CONSTRUIR HISTOGRAMAS
!==========================================================
use globales
use estadistica

implicit none

integer, parameter        :: N=10
real(dp),dimension(N)     :: cuentas,bins
real(dp),dimension(2)     :: lim_op
integer, parameter        :: Npu = 10
real(dp), dimension(Npu)  :: x
integer                   :: i

lim_op = [1_dp,5_dp]

open(20, file ='energia.dat')
do i = 1, Npu
  read(20,*) x(i)
end do
close(20)

call histograma_vec(cuentas,bins,x,N)

write(*,'(2(E12.5))') ( bins(i), cuentas(i), i=1,N)

call histograma(cuentas,bins,'energia.dat',N)

print *, '--------------------------------'
write(*,'(2(E12.5))') ( bins(i), cuentas(i), i=1,N)

end program main 

