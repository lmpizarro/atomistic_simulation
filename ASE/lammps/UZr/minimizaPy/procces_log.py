# -*- coding: utf-8 -*-

import numpy as np

class Process_log(object):

    def __init__ (self, logname):
        self.logname = logname


    def log2numpy(self):
     
        with open (self.logname) as f:
            self.content = f.readlines()
        self.total_lines = len(self.content)

    
        for i,line in enumerate(self.content):
            if line.startswith("Step") == True:
               self.start = i
            if line.startswith("Loop") == True:
               self.end = i
            if line.startswith("Lattice") == True:
               self.lattice = self.content[i] 
            if line.startswith("Cohesive") == True:
               self.cohesive = self.content[i] 
            if line.startswith("Total energy") == True:
               self.energy = self.content[i] 
            if line.startswith("Number of atoms") == True:
               self.natoms = self.content[i] 



        del self.content[self.end-1: self.total_lines]
        del self.content[0:self.start]
 
        print self.cohesive, self.lattice, self.energy, self.natoms

        self.cols = self.content[0].split()
        del self.content[0]
        matrix = []
        
        for l in self.content:
            temp = []
            line = l.split()
            for e in line:
                temp.append(float(e))
            matrix.append(temp)

        self.matrix = np.asarray(matrix)




if __name__ == "__main__":
    pl = Process_log("log")
    matrix = pl.log2numpy()
    print pl.start, pl.end
    print pl.cols
    print pl.matrix[:,1]
