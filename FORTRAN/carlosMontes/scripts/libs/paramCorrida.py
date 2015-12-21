#!/usr/bin/python
# -*- coding: utf-8 -*-


class Corrida(object):
    def __init__(self, Temp, DeltaTiempo, 
            PasosIntegracion, CantidadMediciones, Gamma):

        self.Temp = Temp
        self.Dt = DeltaTiempo
        self.Ntime = PasosIntegracion
        self.Nmed  = CantidadMediciones
        self.Gamma = Gamma

    def write_config(self, name):
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
