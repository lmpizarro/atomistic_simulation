!
!
!
!
!
program main
  use types
  use globales
  use read_param

  implicit none

  call leer_parametros()
  call inicializar_globales()
  call finalizar_globales()

endprogram main
