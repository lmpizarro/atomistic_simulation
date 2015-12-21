#!/usr/bin/python
# -*- coding: utf-8 -*-
from types import *


class Corrida(object):
    def __init__(self, Temp, DeltaTiempo, 
            Nsteps, Nmeasures, Gamma):

        assert (Temp >= 0), "Colder than zero!!"
        assert (DeltaTiempo > 0), "Delta Time mut be greater than 0!!"
        assert (Nsteps > 1000), "Nsteps must be greater than 0!!"
        assert (Nmeasures > 10), "Nmeasures must be greater than 10!!"
        assert (Nmeasures < Nsteps), "Nmeasures must be lower than Nsteps!!"
        assert (Gamma > 0.0), "Gamma must be greater than 0!!"

        assert type(Temp) is FloatType, "Temp must be float"
        assert type(DeltaTiempo) is FloatType, "DeltaTiempo must be float"
        assert type(Nsteps) is IntType, "Nsteps must be integer"
        assert type(Gamma) is FloatType, "Gamma must be float"
        assert type(Nmeasures) is IntType, "Nmeasures must be integer"

        self.Temp = Temp
        self.Dt = DeltaTiempo
        self.Nsteps = Nsteps
        self.Nmed  = Nmeasures
        self.Gamma = Gamma

    def write_config(self, name):

        assert type(name) is StringType, "name must be string"
        end_name = name[-4:]
        assert end_name == ".par", "config file must end with .par"

        format = "%14.8f %14.8f %10d %10d %14.8f\n"
        fo = open(name, 'w')
        tmp = (format)%(self.Temp, self.Dt, self.Ntime, self.Nmed,\
                self.Gamma)
        fo.write(tmp)
        fo.close()

def main ():
    corrida = Corrida(1.0, 0.0001, 10000, 100, 0.5)

if __name__ == "__main__":
    main()
