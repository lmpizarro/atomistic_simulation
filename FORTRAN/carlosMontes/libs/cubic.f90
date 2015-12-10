  !
  ! Para compilar con f2py: f2py -c -m cubic cubic.f90
  !
  subroutine posiciones_fcc(gR, px, py, pz, Dim, Natom )
    implicit none      
 
    integer, intent(in)    :: px, py, pz     ! Periodicidad de la red, # veces que se
                                        ! repite la red en las direcciones x, y e z
    integer :: i, j, k, i_gr, Dim, Natom

    real (kind=8),intent(out), dimension(0:Dim - 1 ,0: Natom - 1 ):: gR
    
!f2py intent(in, out) :: gR  
!f2py depend(Dim) :: gR
!f2py depend(Natom) :: gR
!f2py intent(in) :: px 
!f2py intent(in) :: py 
!f2py intent(in) :: pz 

    i_gr = 0
    do i = 1, px 
      do j= 1,  py
        do k= 1, pz 
           gR(0, i_gr) = DBLE(i - 1)
           gR(1, i_gr) = DBLE(j - 1)
           gR(2, i_gr) = DBLE(k - 1)
           i_gr = i_gr + 1
           gR(0, i_gr) = DBLE(i - 1)
           gR(1, i_gr) = DBLE(j - 0.5)
           gR(2, i_gr) = DBLE(k - 0.5)
           i_gr = i_gr + 1
           gR(0, i_gr) = DBLE(i - 0.5)
           gR(1, i_gr) = DBLE(j - 1)
           gR(2, i_gr) = DBLE(k - 0.5)
           i_gr = i_gr + 1
           gR(0, i_gr) = DBLE(i - 0.5)
           gR(1, i_gr) = DBLE(j - 0.5)
           gR(2, i_gr) = DBLE(k - 1)
           i_gr = i_gr + 1
        enddo
      enddo
   enddo

  end subroutine posiciones_fcc

  subroutine posiciones_bcc(gR, px, py, pz, Dim, Natom )
    implicit none      
 
    integer, intent(in)    :: px, py, pz     ! Periodicidad de la red, # veces que se
                                        ! repite la red en las direcciones x, y e z
    integer :: i, j, k, i_gr, Dim, Natom

    real (kind=8),intent(out), dimension(0:Dim - 1 ,0: Natom - 1 ):: gR
    
!f2py intent(in, out) :: gR  
!f2py depend(Dim) :: gR
!f2py depend(Natom) :: gR
!f2py intent(in) :: px 
!f2py intent(in) :: py 
!f2py intent(in) :: pz 

    i_gr = 0
    do i = 1, px 
      do j= 1,  py
        do k= 1, pz 
           gR(0, i_gr) = DBLE(i - 1)
           gR(1, i_gr) = DBLE(j - 1)
           gR(2, i_gr) = DBLE(k - 1)
           i_gr = i_gr + 1
           gR(0, i_gr) = DBLE(i - 0.5)
           gR(1, i_gr) = DBLE(j - 0.5)
           gR(2, i_gr) = DBLE(k - 0.5)
           i_gr = i_gr + 1
        enddo
      enddo
   enddo

  end subroutine posiciones_bcc

  subroutine posiciones_simple(gR, px, py, pz, Dim, Natom )
    implicit none      
 
    integer, intent(in)    :: px, py, pz     ! Periodicidad de la red, # veces que se
                                        ! repite la red en las direcciones x, y e z
    integer :: i, j, k, i_gr, Dim, Natom

    real (kind=8),intent(out), dimension(0:Dim - 1 ,0: Natom - 1 ):: gR
    
!f2py intent(in, out) :: gR  
!f2py depend(Dim) :: gR
!f2py depend(Natom) :: gR
!f2py intent(in) :: px 
!f2py intent(in) :: py 
!f2py intent(in) :: pz 

    i_gr = 0
    do i = 1, px 
      do j= 1,  py
        do k= 1, pz 
           gR(0, i_gr) = DBLE(i - 1)
           gR(1, i_gr) = DBLE(j - 1)
           gR(2, i_gr) = DBLE(k - 1)
           i_gr = i_gr + 1
        enddo
      enddo
   enddo

  end subroutine posiciones_simple
