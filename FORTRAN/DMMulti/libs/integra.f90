module integra

  use types,             only: dp
  use globales
  use mediciones
  use io_parametros,     only: escribe_trayectoria, escribe_en_columnas

  implicit none

  private

  public       :: integracion, integracion_min, cpc_vec

contains 

  !===============================================================================
  ! integración de las ecuaciones de movimiento - velocity verlet
  !===============================================================================
  ! integra las ecuaciones dinámicas con el algoritmo de velocity-verlet
  subroutine integracion()
    real(dp) :: pres    ! presión instantánea
    real(dp) :: temp    ! presión instantánea
    integer :: Kmed    ! Cantidad de puntos medidos
    real(dp), dimension(:,:), allocatable :: Eng_t   ! energía en función del tiempo
    real(dp), dimension(:), allocatable :: Pres_t  ! presión en función del tiempo
    real(dp), dimension(:), allocatable :: Temp_t  ! temperatura en función del tiempo
    integer :: i, j, k, inic, fin

    ! Se define la cantidad de puntos que se van a medir
    Kmed = int(gNtime/abs(gNmed)) + 1              ! Se agrega +1 para poner el inicial
    allocate( Eng_t(1:3,1:Kmed), Pres_t(1:Kmed), Temp_t(1:Kmed) )

    ! -------------------------------------------------------------------------------------------
    ! COMIENZA EL LOOP PRINCIPAL DE INTEGRACION
    ! ------------------------------------------------------------------------------------------
    k = 1                                              ! Contador para el loop de mediciones
    do i = 1, gNtime 
      inic = 1
      do j=1, gNespecies
        fin = inic + gNp(j) - 1
        !print *, "integracion ", inic, fin, j
        gR(:,inic:fin) = gR(:,inic:fin) + gDt * gV(:,inic:fin) + 0.5_dp * &
                         gF(:,inic:fin) * gDt ** 2 / gLj_param(j,3)   ! gR(t+dt)
        call cpc_vec()
        gV(:,inic:fin) = gV(:,inic:fin) + 0.5_dp * gF(:,inic:fin) * gDt / gLj_param(j,3)     ! gV(t+0.5dt) 
        inic = fin + 1
      enddo

      call calcula_fuerza()                            ! Calcula fuerzas y potencial

      inic = 1
      do j=1, gNespecies
        fin = inic + gNp(j) - 1
        gV(:,inic:fin) = gV(:,inic:fin) + 0.5_dp * gF(:,inic:fin) * gDt / gLj_param(j,3)     ! gV(t+dt)
        inic = fin + 1
      enddo


      ! Se realizan las mediciones
      if (mod(i,gNmed) == 0) then
        ! Energia cinetica
        call calcula_kin()
        ! Temperatura
        call calcula_temp(Temp)
        ! Presión
        call calcula_pres(Pres)
        ! Guarda magnitudes
        Eng_t(:,k) = (/gPot, gKin, gPot + gKin/) 
        Pres_t(k)  = Pres
        Temp_t(k)  = Temp

        print *, "med", k, i, gPot, gKin, gPot + gKin

        k = k + 1                                      ! Actualiza contador mediciones
      end if
    end do
! -------------------------------------------------------------------------------------------
! FIN DEL LOOP PRINCIPAL DE INTEGRACION
! -------------------------------------------------------------------------------------------
    print *, "fin integracion"
    deallocate( Eng_t, Pres_t, Temp_t)
  endsubroutine integracion

  !===============================================================================
  ! INTEGRACIÓN DE LAS ECUACIONES DE MOVIMIENTO - MINIMIZACIÓN ENERGÍA
  !===============================================================================
  ! Subrutina de integración de las ecuaciones de movimiento para minimizar energía
  ! Es el Problema 3 de la Guia_2a
  subroutine integracion_min()

    real(dp), dimension(gNtime+1) :: Eng_t   ! Energía en función del tiempo
    integer :: i, j

    ! Escribe encabezado y primer punto de la trayectoria
    call escribe_trayectoria(gR,.TRUE.) 
    ! El primer punto es la energía inicial
    Eng_t(1) = gPot

    do i = 1, gNtime 
      do j = 1, gNpart
        gR(:,j) = gR(:,j) + 0.5_dp * gF(:,j) * (gDt**2) / gMasa(gIndice_elemento(j))
      end do
      ! Aplica condiciones periódicas de contorno
      call cpc_vec()   

      ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      call escribe_trayectoria(gR,.FALSE.)    
      ! Calcula fuerza y energía
      call calcula_fuerza()
      ! Escribe energía potencial en vector
      Eng_t(i+1) = gPot
    end do

    call escribe_en_columnas(Eng_t/gNpart,'energia_pot_minimizada.dat',gDt)

  end subroutine integracion_min

  subroutine integracion_min_original()

    real(dp), dimension(gNtime+1) :: Eng_t   ! Energía en función del tiempo
    integer :: i, j, inic, fin

    call calcula_fuerza()
    ! Escribe encabezado y primer punto de la trayectoria
    ! call escribe_trayectoria(gR,.FALSE.) 
    ! El primer punto es la energía inicial
    Eng_t(1) = gPot

    do i = 1, gNtime 
      inic = 1
      do j=1, gNespecies
        fin = inic + gNp(j) - 1
        !print *, "integracion ", inic, fin, j
        gR(:,inic:fin) = gR(:,inic:fin) + 0.5_dp * gF(:,inic:fin) * (gDt**2) / gLj_param(j,3)
        inic = fin + 1
      enddo

      ! Aplica condiciones periódicas de contorno
      call cpc_vec()   

      ! Esta subrutine abre y cierra un archivo. Se puede optimizar haciéndolo acá.
      ! call escribe_trayectoria(gR,.FALSE.)    
      ! Calcula fuerza y energía
      call calcula_fuerza()
      ! Escribe energía potencial en vector
      Eng_t(i+1) = gPot
    end do
  end subroutine integracion_min_original

  !===============================================================================
  ! Condiciones períodicas de contorno
  !===============================================================================
  subroutine cpc_vec()

    real(dp), dimension(3,gNpart) :: tmp     ! Variable temporal
    integer :: j
    
    !gR = gR - gLado_caja*floor(gR/gLado_caja)
    
    ! Lo escribo de esta forma porque de lo contrario da error de compilación
    ! en el cluster (dice ser un bug de gfortran)
    tmp = gLado_caja*floor(gR/gLado_caja)

    gR = gR - tmp

  end subroutine cpc_vec

end module integra
