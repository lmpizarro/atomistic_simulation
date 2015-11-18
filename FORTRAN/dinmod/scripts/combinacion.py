
# coding: utf-8

# In[16]:

import numpy as np

# Dimensión del problema
Rho = .1
L = 5 # veces sigma_a
percent_a = 0.3
percent_b = 0.5

Npart = np.int(Rho  * L * L * L)
Na = int(percent_a * Npart)
Nb = int(percent_b * Npart)
Nc = Npart - Na - Nb

Np = np.asarray([Na,Nb, Nc])

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

Epsilon = calc_matriz_epsilon(epsilon) / epsilon[0]
Sigma =  calc_matriz_sigma(sigma) / sigma[0]

print Epsilon
print Sigma


Ma = 1; Eps_a = 1.83592*10E-2; sigma_a = 4.1
Mb = .8; Eps_b = 1.03413*10E-2; sigma_b = 4.0
Mc = .7; Eps_c = 1.2E-2; sigma_c = 3.9

# parámetros de la interacción de átomos diferentes
Eps_ab = np.sqrt(Eps_a * Eps_b); sigma_ab = (sigma_a + sigma_b) / 2
Eps_ac = np.sqrt(Eps_a * Eps_c); sigma_ac = (sigma_a + sigma_c) / 2
Eps_bc = np.sqrt(Eps_b * Eps_c); sigma_bc = (sigma_b + sigma_c) / 2

# Transformación a por unidad, se toma como referencia el átomo A
 
P = np.asarray([i for i in range((Na + Nb) * 3)]).reshape(Na + Nb, 3)

#inicializa posiciones
P = np.random.rand((Na + Nb)*3).reshape(Na + Nb, 3) * L
V = np.random.rand((Na + Nb)*3).reshape(Na + Nb, 3) - 0.5
# el calculo de la velocidad del centro de masa debe tener 
# en cuenta la masa diferente de cada átomo
# vcm = (Ma * Sum (vai) + Mb * Sum (vbi) ) / (Na * Ma + Nb * Mb)
Vcm = (Ma * np.sum(V[:Na], 0) + Mb * np.sum(V[Na:], 0)) /  (Na * Ma + Nb * Mb)


#print (V)
#print (np.sum(V[:Na], 0))
#print (np.sum(V[Na:], 0))
print (Vcm)


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



