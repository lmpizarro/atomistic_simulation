#!/usr/bin/python
# -*- coding: utf-8 -*-



import paramCorrida as PC
import paramMezcla as PM
import cubic as CB


def main():
    # inicializo los componentes del cálculo
    components = PM.Interacciones_LJ(1.2, 4.0, 2.3, "componente1")
    components.add_component(1.1, 4.5, 2.1, "componente2")
    components.add_component(1.0, 4.3, 2.0, "componente2")
    components.calculate_matrices()
    # inicializo la estructura
    fcc = CB.Cubic(2, nx=4, ny=4, nz=4, tipo="fcc")
    fcc.set_A2B2()
    # parámetros de la corrida
    corrida = PC.Corrida(1.0, 0.0001, 10000, 100, 0.5)
    name_config = "config.par"
    corrida.write_config(name_config)
    components.write_file(name_config)
    fcc.write_config(name_config)


if __name__ == "__main__":
    main()
