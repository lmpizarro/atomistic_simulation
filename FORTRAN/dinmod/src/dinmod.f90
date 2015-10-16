program main_dimod 

    use globales,         only: gT, gDt, gL, gNpart, gNtime, gSigma, gEpsil, gM
    use dinmods,          only: inicializacion, finalizacion, integracion_min, integracion 

    implicit none

    ! Inicializa parámetros 
    call inicializacion()
    ! Busca el mínimo de energía
    call integracion_min()
    ! Integración de la dinámica
    call integracion()
    ! Finalización del programa
    call finalizacion()

end program main_dimod
