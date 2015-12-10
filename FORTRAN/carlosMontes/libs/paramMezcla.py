#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

class Interacciones_LJ(object):
    def __init__(self, epsilon, sigma, masa, name):
        self.Componentes = [{"name": name, "parameters": [epsilon, sigma,
            masa]}]

    def add_component(self, epsilon, sigma, masa, name):
        c = {"name": name, "parameters": [epsilon, sigma, masa]}
        self.Componentes.append(c)
        self.Ncomponents = len(self.Componentes)

    def list_components(self):
        for i, c in enumerate(self.Componentes):
            print "Elemento: ", i + 1
            print "Nombre: ", c["name"]
            print "Epsilon: ", c["parameters"][0]
            print "Sigma: ", c["parameters"][1]
            print "Masa: ", c["parameters"][2]
            print "-----------------"

    def calculate(self):
        self.sigma = np.zeros(len(self.Componentes) **
                2).reshape(self.Ncomponents, self.Ncomponents)

        self.epsilon = np.zeros(len(self.Componentes) **
                2).reshape(len(self.Componentes), len(self.Componentes))


        for i, c in enumerate(self.Componentes):
            self.epsilon[i][i] =  c["parameters"][0]
            self.sigma[i][i] =  c["parameters"][1]


        for i in range(self.Ncomponents):
            for j in range(i + 1, self.Ncomponents):
                self.sigma[i][j] = (self.Componentes[i]["parameters"][1] + 
                                self.Componentes[j]["parameters"][1]) /2
                self.sigma[j][i] = self.sigma[i][j]

        #normalizacion
        self.sigma = self.sigma / self.sigma[0][0]

        for i in range(self.Ncomponents):
            for j in range(i + 1, self.Ncomponents):
                self.epsilon[i][j] = np.sqrt((self.Componentes[i]["parameters"][0] * 
                                self.Componentes[j]["parameters"][0]))
                self.epsilon[j][i] = self.epsilon[i][j]

        #normalizacion
        self.epsilon = self.epsilon / self.epsilon[0][0]

    def list_matrices(self):
        print "Sigma"
        print ""
        print self.sigma
        print ""
        print "Epsilon"
        print self.epsilon
        print ""

def main ():
    components = Interacciones_LJ(1.2, 4.0, 2.3, "componente1")
    components.add_component(1.1, 4.5, 2.1, "componente2")
    components.add_component(1.0, 4.3, 2.0, "componente2")
    components.add_component(1.05, 4.4, 2.7, "componente2")
    components.list_components()
    components.calculate()
    components.list_matrices()

if __name__ == "__main__":
    main()
