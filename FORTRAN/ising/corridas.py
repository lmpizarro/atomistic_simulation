#!/usr/bin/env python


tempe = [0.5, 0.6, 0.7, 1.2, 1.4]
tempe = [str(i) for i in tempe]

for T in tempe:
    with open('parametros.dat','r') as f1:
        for line in f1:
            col = line.split()
            new=line.replace(col[2],T)
   
    with open('parametros.dat','w') as fileout:
        fileout.write(new)

    # Creo directorio
    # copio ejecutable y entradas
    # Hago algo con la salida
    