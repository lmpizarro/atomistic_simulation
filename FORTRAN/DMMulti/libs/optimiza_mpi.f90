!
! modulo que tiene como fin 
! crear funciones para balancear
! carga entre procesadores
!
module optimiza_mpi

  implicit none


contains
  ! recibe un array de tamaño size + 1, donde
  ! size es la cantidad de procesadores 
  ! recibe la cantidad de partículas sobre la
  ! que se hace la DM
  ! para 256 o mas particulas
  subroutine optimiza_para(pasos_r, gNpart)
    integer, intent(inout), dimension(:) :: pasos_r
    integer :: size_, i, j, gNpart
    real :: area, tmp
    
    size_ = size(pasos_r) ! la cantidad de procesadores + 1
    area =  (gNpart * gNpart) / (2*(size_ - 1))

    pasos_r(1) = 1
    do j = 1, size_ - 1
      do i = 1, gNpart
        tmp = (i - pasos_r(j))*((i -pasos_r(j)) / 2 - i + gNpart)
        if (tmp .ge. area) then
          pasos_r(j + 1) = i
          print *, "mmmmmmmmmm", i
          !inic = i
          exit
        endif   
      enddo
    enddo

    pasos_r(size_ ) = gNpart

  endsubroutine optimiza_para

  ! call test()
  subroutine test()
    integer :: gNpart
    integer :: inic, fin

    integer :: i, j 
    integer :: size_, rank, k

    integer, allocatable:: limites_procesador (:)

    integer, allocatable:: cantidadCuentas (:)

    size_ = 10

    gNpart = 500

    allocate (limites_procesador(1:size_ + 1))


    call optimiza_para(limites_procesador, gNpart)

    print *, "despues de optimiza"
    do i=1, size_  + 1
      print *, limites_procesador(i), size(limites_procesador)
    enddo


    fin= 0
    inic = limites_procesador(1)
    do k=1, size_
      fin = limites_procesador(k + 1)
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
  endsubroutine test
          
endmodule optimiza_mpi
