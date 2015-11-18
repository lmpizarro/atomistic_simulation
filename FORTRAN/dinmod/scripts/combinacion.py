
# coding: utf-8

# In[46]:

import numpy as np


#
# Calcula la cantidad de partículas de un sistema de Nespecies
# y la cantidad de partículas por especies dados los porcentajes
# atómicos de Nespecies - 1
#
Rho = .9
L = 5 # veces sigma_a
percent_a = 0.3
percent_b = 0.5

Npart = np.int(Rho  * L * L * L)
Na = int(percent_a * Npart)
Nb = int(percent_b * Npart)
Nc = Npart - Na - Nb

Np = np.asarray([Na,Nb, Nc])

print "cantidad de partículas total %d" % (Npart)

for i in range(np.size(Np)):
    print "cantidad de partículas del elemento %d: %d" % (i, Np[i])

    
#
# Calcula los parámetros de interacción de un sistema
# de Nespecies, dados los vectores epsilon y sigma
#    
#epsilon en eV
epsilon = np.asarray([1.83592*10E-2, 1.03413*10E-2, 1.2*10E-2])
#sigma en angstrom
sigma   = np.asarray([4.1,4.0,3.9])
masa = np.asarray([1,.8,.7])
Nespecies = np.size(epsilon)

def calc_matriz_epsilon(vec_parametros):
    Nespecies = np.size(vec_parametros)
    matriz_interaccion = np.zeros(Nespecies*Nespecies).reshape(Nespecies,Nespecies)
    
    for i in range(Nespecies):
        matriz_interaccion[i,i] = vec_parametros[i]

    for i in range(np.size(vec_parametros) - 1):
        for j in range (i + 1,Nespecies):
            matriz_interaccion[i,j] = np.sqrt(vec_parametros[i] * vec_parametros[j])
            matriz_interaccion[j,i] =  matriz_interaccion[i,j]
    return matriz_interaccion

def calc_matriz_sigma(vec_parametros):
    Nespecies = np.size(vec_parametros)
    matriz_interaccion = np.zeros(Nespecies*Nespecies).reshape(Nespecies,Nespecies)
    
    for i in range(Nespecies):
        matriz_interaccion[i,i] = vec_parametros[i]

    for i in range(np.size(vec_parametros) - 1):
        for j in range (i + 1,Nespecies):
            matriz_interaccion[i,j] = (vec_parametros[i]  + vec_parametros[j])/2
            matriz_interaccion[j,i] =  matriz_interaccion[i,j]
    return matriz_interaccion


# Cálculo y Transformación a por unidad, 
# se toma como referencia el átomo de la posición 0
Epsilon = calc_matriz_epsilon(epsilon) / epsilon[0]
Sigma =  calc_matriz_sigma(sigma) / sigma[0]


#inicializa posiciones
P = np.random.rand((Npart)*3).reshape(Npart, 3) * L
V = np.random.rand((Npart)*3).reshape(Npart, 3) - 0.5
# el calculo de la velocidad del centro de masa debe tener 
# en cuenta la masa diferente de cada átomo
# vcm = (Ma * Sum (vai) + Mb * Sum (vbi) ) / (Na * Ma + Nb * Mb)


Vcm = np.zeros(3)
inic = 0; fin = 0
for i in range(np.size(Np)):
    fin = inic + Np[i]
    Vcm = Vcm + masa[i] * np.sum(V[inic:fin],0)
    inic = Np[i]

Vcm = Vcm / np.sum(Np * masa)
print "Velocidad del centro de masa: %.3f %.3f %.3f" % (Vcm [0], Vcm[1], Vcm[2])


# In[58]:

def calc_energ(P, Np):
    
    cant_atomos = np.size(Np)
 
    #print P
    
    print "calcula interacciones entre iguales"
    inic = 0
    for i in range(cant_atomos):
        fin = inic +  Np[i]
        print inic, fin 
        #print P[inic:fin]
        inic = fin
        


# In[59]:

calc_energ(P, Np)


# In[ ]:



