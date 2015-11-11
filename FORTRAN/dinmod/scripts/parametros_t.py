# -*- coding: utf-8 -*-

###############################################################################       
#   PARAMETROS DE ENTRADA PARA CORRER A DISTINTAS TEMPERATURAS
###############################################################################
# Número de partículas
N_part = 200
# Densidad de partículas
Rho = 0.3
#------ Barrido de temperaturas
# Temperatura mínima
Temp_min = 0.7
# Temperatura máxima
Temp_max = 1.4
# Paso de temperatura
dTemp = 0.25
# Agrego el detalle cerca de la temperatura crítica
#T_detail_min = 2.10
#T_detail_max = 2.50
#dT_detail = 0.02
# abs(K_grab) Cada cuántos puntos se quiere grabar el archivo temporal
# K_grab < 0 especifica que no se grabe ningún archivo temporal
N_grab = 10
# Paso temporal de integración
dt = 0.001
# Número de pasos para la primer corrida (termalización)
N_term = 1000
# Número de pasos para la segunda corrida (medición)
N_medi = 2000
# Número de corridas para cada temperatura
Nrun = 8
#
# FIN PARAMETROS DE ENTRADA
###############################################################################


