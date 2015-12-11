program md_read
  use param_md, only: Parametros

  implicit none


  type(Parametros) :: params

  call params % leer()
endprogram
