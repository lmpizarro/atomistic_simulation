# -*- coding: utf-8 -*-

###############################################################################       
#   PARAMETROS DE ENTRADA PARA CORRER A DISTINTAS TEMPERATURAS
###############################################################################
# Número de partículas
N_part = 200
# Temperatura
Temp = 1.1
# Densidad de partículas
Rho = [0.001, 1.0]
#Rho = [0.001, 0.01, 0.1, 0.8, 0.9, 1.0]
#------ Barrido de densidades
# Densidad mínima
#Rho_min = 0.7
# Densidad máxima
#Rho_max = 1.4
# Paso de Densidad
#dRho = 0.25
# Agrego el detalle 
Rho_detail_min = 0.10
Rho_detail_max = 0.80
dRho_detail = 0.05
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


