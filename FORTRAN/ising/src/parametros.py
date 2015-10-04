# -*- coding: utf-8 -*-

###############################################################################       
#   PARAMETROS DE ENTRADA
###############################################################################
# Tamaño de la red de spines
N_red = 20
M_red = 20
# Cada cuántos puntos se quiere grabar el archivo temporal
# K_grab = 0 especifica que no se grabe ningún archivo temporal
N_grab = 0     
# Barrido de temperaturas
# Temperatura mínima
T_min = 0.5
# Temperatura máxima
T_max = 4.0
# Paso de temperatura
dT = 0.1
# Agrego el detalle cerca de la temperatura crítica
T_detail_min = 2.15
T_detail_max = 2.40
dT_detail = 0.02
# Número de pasos para la primer corrida (termalización)
N_term = '4000'
# Número de pasos para la segunda corrida (medición)
N_medi = '10000'
# Número de corridas para cada temperatura
Nrun = 8
#
# FIN PARAMETROS DE ENTRADA
###############################################################################


