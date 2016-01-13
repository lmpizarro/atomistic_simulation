# -*- coding: utf-8 -*-

class MEAMPOT_UZr(object):
    potentials =[{"type":"meam", "name": "UZr", "params":{"meafile" :'''rc = 5.5
delta(1,2) = 0.8
re(1,2) = 2.85
delr = 0.1
alpha(1,2) = 4.800000
lattce(1,2) = b1
rho0(1) = 1.2 
rho0(2) = 0.8
ialloy = 1
emb_lin_neg = 1
bkgd_dyn = 1
zbl(1,1) = 1000
zbl(1,2) = 1000
zbl(2,2) =1000
attrac(1,1) = 0.105
repuls(1,1) = 0.105
attrac(1,2) = 0.05
repuls(1,2) = 0.05
attrac(2,2) = -0.03
repuls(2,2) = -0.03
Cmin(1,1,1) = 1.00
Cmax(1,1,1) = 1.9
Cmin(1,1,2) = 0.5
Cmax(1,1,2) = 2.0
Cmin(1,2,1) = 0.8
Cmax(1,2,1) = 2.5
Cmin(1,2,2) = 0.5
Cmax(1,2,2) = 2.8
Cmin(2,2,1) = 0.5
Cmax(2,2,1) = 2.8
Cmin(2,2,2) = 0.7 
Cmax(2,2,2) = 0.99''',

    "meamf":'''# meam data from vax files fcc,bcc,dia
elt        lat     z       ielement     atwt
alpha      b0      b1      b2           b3    alat    esub    asub
t0         t1              t2           t3            rozero  ibar
'U'        'fcc'    12      92        238.03             
5.1         4.80    6.0     6         6        4.280   5.27   0.98
1           2.50            4         1.0       1       0   
'Zr'       'bcc'    8.      40        91.224
4.10        2.80    2.0     7.00      1        3.535   6.20   0.48
1.0         3.00            2.0      -7       1       0'''}}, 
{"type":"meam", "name": "UNb", "params":{}}, 
{"type":"eam", "name": "UO", "params":{}}]

    def __init__(self, type, name):
        self.type = type
        self.name = name
        self.params = None 
        self.index = None

        for i, p in enumerate(self.potentials):
            if p["type"] == self.type:
                if p["name"] == self.name:
                    self.index = i

        if self.index != None:
            self.params = self.potentials[self.index]  

    def gen_files(self):
        if pot.params != None:
            print pot.params


        fout = open('meafile','w')
        fout.write(self.meafile)
        fout.close()

        fout = open('meamf','w')
        fout.write(self.meamf)
        fout.close()

if __name__ == "__main__":
    # graba los potenciales de acuerdo a como los espera el infile
    pot =  MEAMPOT_UZr("meam", "UZr")

