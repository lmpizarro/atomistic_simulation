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
Temp_min = 200 
# Temperatura máxima
Temp_max = 2000 
# Paso de temperatura
dTemp =  50 
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
N_term = '2000'
# Número de pasos para la segunda corrida (medición)
N_medi = '5000'
# Número de corridas para cada temperatura
Nrun = 1 
epsilon = 1.0
sigma   = 1.0
masa = 1.0
N_save  = 1000
gamma   = 0.5

Size = (3,3,3)
Pbc  = True
LatticeConstant = 3.5

#
# FIN PARAMETROS DE ENTRADA
###############################################################################


