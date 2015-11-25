program main

  implicit none

  integer :: gNpart, inic, fin
  integer :: i, j, size_, rank, k
  integer, allocatable:: delta_x (:)
  real :: coef1, coef2
  real :: num, sum_
  integer, allocatable:: cantidadCuentas (:)

  size_ = 10
  coef2 = 0.9
  coef1 = 2.0
  gNpart = 256

  allocate (delta_x(1:size_))

  do i=1, size_
     num = i
     delta_x(i) =  (num ** coef2 ) + coef1
  enddo

  sum_ = sum(delta_x)

  print *, "suma ", sum_

  delta_x = int(gNpart * delta_x / sum_)

  do i = 1, size_
    print *, "delta x", delta_x(i)
  enddo

  fin= 0
  inic = 1
  do k=1, size_
    fin = fin + delta_x(k)
    if (k.eq.size_) then
      if (fin.lt.gNpart) then
        fin = gNpart
      endif
    endif
    !print *, inic, fin, int(delta_x(k))
     
    rank = 0
    do i=inic, fin - 1
      do j = i + 1 , gNpart 
        rank = rank + 1
        !print *, i, j
      enddo
    enddo
    print * , "total del cuentas: ", rank
    inic = fin
  enddo   
endprogram main
