#!/usr/bin/python
# -*- coding: utf-8 -*-

from types import *
import numpy as np


#
#  clase usada para el cálculo de las matrices de interacción
#
class Interacciones_LJ(object):
    # inicializa la clase con un componente


    def __init__(self, epsilon, sigma, masa, name):
        self.assert_inputs_component(epsilon, sigma, masa, name)
        self.Componentes = [{"name": name, "parameters": [epsilon, sigma,
                                                          masa]}]

    def assert_inputs_component(self, epsilon, sigma, masa, name):
        assert type(epsilon) is FloatType, "epsilon must be float"
        assert type(sigma) is FloatType, "sigma must be float"
        assert type(masa) is FloatType, "masa must be float"
        assert type(name) is StringType, "masa must be float"

        assert (epsilon > 0), "epsilon mut be > 0"
        assert (sigma > 0), "sigma mut be > 0"
        assert (masa > 0), "masa mut be > 0"


    # agrega componentes
    def add_component(self, epsilon, sigma, masa, name):
        self.assert_inputs_component(epsilon, sigma, masa, name)

        c = {"name": name, "parameters": [epsilon, sigma, masa]}
        self.Componentes.append(c)
        self.Ncomponents = len(self.Componentes)

    # lista los componentes de la estructura
    def list_components(self):
        for i, c in enumerate(self.Componentes):
            print "Elemento: ", i + 1
            print "Nombre: ", c["name"]
            print "Epsilon: ", c["parameters"][0]
            print "Sigma: ", c["parameters"][1]
            print "Masa: ", c["parameters"][2]
            print "-----------------"

    # calcula las matrices
    def calculate_matrices(self):
        self.sigma = np.zeros(len(self.Componentes) **
                              2).reshape(self.Ncomponents, self.Ncomponents)

        self.epsilon = np.zeros(
            len(self.Componentes) ** 2).reshape(len(self.Componentes), len(self.Componentes))

        self.masa = np.zeros(len(self.Componentes))
        self.potCut = np.zeros(
            len(self.Componentes) ** 2).reshape(len(self.Componentes), len(self.Componentes))

        for i, c in enumerate(self.Componentes):
            self.epsilon[i][i] = c["parameters"][0]
            self.sigma[i][i] = c["parameters"][1]
            self.masa[i] = c["parameters"][2]

        self.masa = self.masa / self.masa[0]

        for i in range(self.Ncomponents):
            for j in range(i + 1, self.Ncomponents):
                self.sigma[i][j] = (self.Componentes[i]["parameters"][1] +
                                    self.Componentes[j]["parameters"][1]) / 2
                self.sigma[j][i] = self.sigma[i][j]

        # normalizacion
        self.sigma = self.sigma / self.sigma[0][0]

        for i in range(self.Ncomponents):
            for j in range(i + 1, self.Ncomponents):
                self.epsilon[i][j] = np.sqrt(
                    (self.Componentes[i]["parameters"][0] *
                     self.Componentes[j]["parameters"][0]))
                self.epsilon[j][i] = self.epsilon[i][j]

        # normalizacion
        self.epsilon = self.epsilon / self.epsilon[0][0]

        self.calc_radio_corte()
        self.calc_pote_despla()

    # falta calcular el radio de corte y el desplazamiento

    # calcula el radio de corte
    def calc_radio_corte(self):
        self.rc2 = (2.5 * self.sigma) ** 2

    # calcula el potencial de desplazamiento
    def calc_pote_despla(self):
        for i in range(self.Ncomponents):
            self.potCut[i][i] = 4.0 * self.epsilon[i][i] * \
                (self.rc2[i][i] ** (-6) - self.rc2[i][i] ** (-3))

        for i in range(self.Ncomponents):
            for j in range(self.Ncomponents):
                self.potCut[i][j] = 4.0 * self.epsilon[i][j] * \
                    (self.rc2[i][j] ** (-6) - self.rc2[i][j] ** (-3))
                self.potCut[j][i] = self.potCut[i][j]

    # lista las matríces
    def list_matrices(self):
        print "Sigma"
        print ""
        print self.sigma
        print ""
        print "Epsilon"
        print self.epsilon
        print ""
        print "Masa"
        print self.masa
        print ""
        print "RC2"
        print self.rc2
        print ""
        print "PotCut"
        print self.potCut

    def __repr__(self):
	sig = "sigma:\n"
        format = "%14.8f"
        for s in self.sigma:
            for n in s:
                sig = sig + (format) % (n) + "  "
            sig = sig + "\n"

	eps = "epsilon:\n"
        for s in self.epsilon:
            for n in s:
                eps = eps + (format) % (n) + "  "
            eps = eps + "\n"

	rc = "rc2: \n"
        for s in self.rc2:
            for n in s:
                rc = rc + (format) % (n) + "  "
            rc = rc + "\n"

	pc = "potcut: \n"
        for s in self.potCut:
            for n in s:
                pc = pc + (format) % (n) + "  "
            pc = pc + "\n"

	ma = "masses: \n"
        for n in self.masa:
            ma = ma + (format) % (n) + "  "
            ma = ma + "\n"

        total = ("Nelements %4d\n") % (self.Ncomponents) + eps + sig + rc + pc + ma

        return total

    def write_file(self, name):

        assert type(name) is StringType, "name must be string"
        end_name = name[-4:]
        assert end_name == ".par", "config file must end with .par"

        f = open(name, 'a')
        f.write(self.toString())
        f.close

# tests
def main():
    components = Interacciones_LJ(1.2, 4.0, 2.3, "componente1")
    components.add_component(1.1, 4.5, 2.1, "componente2")
    components.add_component(1.0, 4.3, 2.0, "componente2")
    components.add_component(1.05, 4.4, 2.7, "componente2")
    components.list_components()
    components.calculate_matrices()
    components.list_matrices()
    print components.toString()
    components.write_file("config.par")

if __name__ == "__main__":
    main()
