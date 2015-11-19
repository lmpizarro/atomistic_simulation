!
!
!
!
!
program main
  use types
  use globales
  use read_param
  use combinacion

  implicit none

  call leer_parametros()
  call inicializar_globales()
  call comb_sigma()
  call comb_epsilon()
  call finalizar_globales()

endprogram main
