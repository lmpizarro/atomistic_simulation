#!/usr/bin/env python
from numpy import *
import argparse,sys

#
# uso para lammps
# $ python create_crystal.py -s 2 -a 4.5 -n 16 16 16 -t U -d Zr -o 2
#


examples = "EXAMPLES:\n\
create_crystal.py -s1   -a 1.50       -n 4 4 4  # simple cubic,\n\
create_crystal.py -s2   -a 2.20       -n 3 3 3  # BCC,\n\
create_crystal.py -s3   -a 2.20       -n 4 4 4  # FCC / rock-salt,\n\
create_crystal.py -s4   -a 1.42 3.30  -n 4 4 2  # graphite,\n\
create_crystal.py -s4   -a 1.42       -n 4 4 1  # graphene,\n\
create_crystal.py -s5   -a 3.46       -n 2 2 2  # diamond / zinc-blende,\n\
create_crystal.py -s6   -a 2.50 4.00  -n 4 4 4  # wurtzite"
parser = argparse.ArgumentParser(description="create a crystal structure file",epilog=examples,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-s", dest="structure",type=int, default=1, help="choose structure (1=sc, 2=bcc, 3=fcc, 4=graphite, 5=diamond, 6=wurtzite)")
parser.add_argument('-n', metavar='i', dest="ncells", type=int, nargs=3, default = [1, 1, 1],  help="number of unit cells in x, y and z direction (e.g. -n 3 2 4)")
parser.add_argument('-a', metavar='x', dest="lattice_constants", type=float, nargs='+', default = [1], help="supply the necessary lattice constants")
parser.add_argument('-t', dest="atomtype",type=str, default="C", help="primary atom type")
parser.add_argument('-d', dest="diatomic",type=str, help="second atom type (creates diatomic structure)")
parser.add_argument('-b', dest="basename", default="base", help="basename for output file")
parser.add_argument("-o", dest="output",type=int, default=1, help="type of output (1=xyz,2=lammps,3=for,4=(x)bs),5=obj")

def main():
    args      = parser.parse_args()
    lc        = args.lattice_constants
    basename  = args.basename

    if args.structure is 1:
      print "SIMPLE CUBIC"; a=lc[0]
      cell,coords = sc(a)
    elif args.structure is 2:
      print "BCC"; a=lc[0]
      cell,coords = bcc(a)
    elif args.structure is 3:
      print "FCC"; a=lc[0]
      cell,coords = fcc(a)
    elif args.structure is 4 and len(lc)==1:
        print "GRAPHENE"; rcc=lc[0]; 
        cell,coords = graphene(rcc)
        # expand to orthorhombic cell
        cell,coords = expand((1,2,1),cell,coords)
        cell=(array([2*cell[0][0],0,0]),array([0,cell[1][1],0]),array([0,0,cell[2][2]]))
    elif args.structure is 4 and len(lc)==2:
        print "GRAPHITE"; rcc = lc[0]; d = lc[1]
        cell,coords = graphite(rcc,d)
        # expand to orthorhombic cell
        cell,coords = expand((1,2,1),cell,coords)
        cell=(array([2*cell[0][0],0,0]),array([0,cell[1][1],0]),array([0,0,cell[2][2]]))
    elif args.structure is 5:
      print "DIAMOND"; a=lc[0]
      cell,coords = diamond(a,array([0,0,0]))
    elif args.structure is 6 and len(lc)==2:
      print "WURTZITE"; a=lc[0]; c=lc[1]
      cell,coords = wurtzite(a,c)
    else:
      print 'structure not implemented or incorrect number of lattice constants'; return
  
    # expand crystal (n,m,k) times
    cell,coords = expand(args.ncells,cell,coords)

    # define atom types
    natoms = len(coords)
    types = [args.atomtype]*natoms
    if args.diatomic: 
      for i in range(1,natoms,2): types[i] = args.diatomic

    # compute PBCs
    (a1,a2,a3) = cell
    a = sqrt(dot(a1,a1)); 
    b = sqrt(dot(a2,a2)); 
    c = sqrt(dot(a3,a3)); 
    if(a==0): a=1
    if(b==0): b=1
    if(c==0): c=1
    print a,b,c
    alpha = arccos(dot(a2,a3)/(b*c))*180/pi
    beta  = arccos(dot(a1,a3)/(a*c))*180/pi
    gamma = arccos(dot(a1,a2)/(a*b))*180/pi
    print dot(a1,a2),(a*b)
    print dot(a1,a2)/(a*b)
    print arccos(dot(a1,a2)/(a*b))
    print arccos(dot(a1,a2)/(a*b))*180/pi
    print '-----------------------------------------'
    print 'created ',len(coords),'atom sample'
    print '(super)cell vectors'
    print "a1 = array([% 6.2f, % 6.2f, % 6.2f])" % (a1[0],a1[1],a1[2])
    print "a2 = array([% 6.2f, % 6.2f, % 6.2f])" % (a2[0],a2[1],a2[2])
    print "a3 = array([% 6.2f, % 6.2f, % 6.2f])" % (a3[0],a3[1],a3[2])
    print '-----------------------------------------'
    print "cell for VMD : pbc set { %.3f %.3f %.3f   %.2f %.2f %.2f }; pbc box"%(a,b,c,alpha,beta,gamma)

    if   args.output is 1: writexyz(basename+".xyz",cell,coords,types)
    elif args.output is 2: writelammps(basename+".data",cell,coords, types)
    elif args.output is 3: writefor(basename+".for",basename+".sys",cell,coords,types)
    elif args.output is 4: writebs(basename+".bs",cell,coords,types)
    elif args.output is 5: writeobj(basename+".obj",coords)
    else: print 'output method not implemented'; return


######### crystal structure functions ##########
def sc(a,r=array([0,0,0])): # simple cubic
    r0 = r
    cell = a*identity(3)
    return cell,[r0]

def bcc(a,r=array([0,0,0])): # body centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,1])
    cell = a*identity(3)
    return cell,[r0,r1]

def fcc(a,r=array([0,0,0])): # face centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,0])
    r2 = r + dot((a/2.),[1,0,1])
    r3 = r + dot((a/2.),[0,1,1])
    cell = a*identity(3)
    return cell,[r0,r1,r2,r3]

def graphene(rcc,r=array([0,0,0])): # rcc = CC distance in graphene
    # translation vectors for generating graphene:
    a1 = dot(rcc,[3/2.,-sqrt(3)/2,0])
    a2 = dot(rcc,[3/2.,+sqrt(3)/2,0])
    a3 = array([0,0,0])
    cell = (a1,a2,a3)

    # two atoms per unit cell
    t = array([rcc,0,0])
    r0 = r 
    r1 = r + t
    coords = [r0,r1]
    return cell,coords

def graphite(rcc,d,r=array([0,0,0])): # rcc = CC distance in graphene, d = interplanar distance
    # translation vectors for generating graphite:
    a1 = dot(rcc,[3/2.,-sqrt(3)/2,0])
    a2 = dot(rcc,[3/2.,+sqrt(3)/2,0])
    a3 = dot(d,[0,0,2])
    cell = (a1,a2,a3)

    # four atoms per unit cell
    t0 = dot(rcc,[1,0,0])
    t1 = dot(d,[0,0,1])
    r0 = r
    r1 = r0 + t0
    r2 = r1 + t1
    r3 = r2 + t0
    coords = [r0,r1,r2,r3] 
    return cell,coords

def diamond(a,r=array([0,0,0])):
    # diamond (or zinc-blende, = diamond with alternating atoms)   r_CC = sqrt(3)/4 * a        (3.46)
    r0 = r 
    r1 = r + dot((a/4.),[1,1,1])
    r2 = r + dot((a/4.),[2,2,0])
    r3 = r + dot((a/4.),[3,3,1])
    r4 = r + dot((a/4.),[0,2,2])
    r5 = r + dot((a/4.),[3,1,3])
    r6 = r + dot((a/4.),[2,0,2])
    r7 = r + dot((a/4.),[1,3,3])
    coords = [r0,r1,r2,r3,r4,r5,r6,r7]
    cell = a*identity(3)
    return cell,coords

def wurtzite(a,c,r=array([0,0,0])):
    # from http://cst-www.nrl.navy.mil/lattice/struk/b4.html
    u=3./8. 

    # primitive vectors
    a1 = dot(a/2.,[1,-sqrt(3),0])
    a2 = dot(a/2.,[1, sqrt(3),0])
    a3 = dot(c   ,[0,0,1])
    cell = (a1,a2,a3)

    # basis vectors
    b1 = dot(1/3.,a1) + dot(2/3.,a2)
    b2 = dot(1/3.,a1) + dot(2/3.,a2) + dot(u,a3)
    b3 = dot(2/3.,a1) + dot(1/3.,a2) + dot(1/2.,a3)
    b4 = dot(2/3.,a1) + dot(1/3.,a2) + dot(1/2.+u,a3)

    coords = [b1,b2,b3,b4]
    return cell,coords

######## Misc. Functions #########
def expand((n1,n2,n3),cell,coords): 
    # expand current cell and coordinates by a1,a2,a3
    newcoords = []
    (a1,a2,a3) = cell
    for x in range(n1):
        for y in range(n2):
            for z in range(n3):
                r = dot(x,a1) + dot(y,a2) + dot(z,a3)
                for c in coords:
                    newcoords.append(c + r)

    newcell = (dot(n1,a1),dot(n2,a2),dot(n3,a3))
    return newcell,newcoords

def periodic(cell,coords):
    # move atoms into periodic cell box
    for c in coords:
        for i in range(3):
            c[i] = c[i] % cell[i]

def translate(r,coords):
    coords = [r+c for c in coords]

# Input/Output Functions
def writexyz(file,cell,coords,types):
    # xyz
    f=open(file,'w')
    f.write(str(len(coords))+'\n')
    f.write("%12.8f %12.8f %12.8f\n" % (cell[0][0], cell[1][1], cell[2][2]))
    for i,l in enumerate(coords):
        f.write("%-2s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))  # xyz format

def writelammps(file,cell,coords, types): 
    # LAMMPS
    elements = unique(types)
    N  = len(elements) 
    f=open(file,'w')
    f.write('info: \n')
    f.write('\n')
    f.write(str(len(coords))+' atoms\n')
    f.write('%d atom types\n'%(N))
    f.write('\n')
    f.write("0.0 %12.8f xlo xhi" % float(cell[0][0]) + '\n')
    f.write("0.0 %12.8f ylo yhi" % cell[1][1] + '\n')
    f.write("0.0 %12.8f zlo zhi" % cell[2][2] + '\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for i,l in enumerate(coords):
        tipo = i % N  + 1
        f.write("%3d %3d %12.8f %12.8f %12.8f\n" % (i+1, tipo, l[0], l[1], l[2])) # lammps format
    f.close()

def writefor(file,sysfile,cell,coords,types):
    # tersoff
    b=0.5291772083
    g=open(sysfile,'w')
    g.write("&SYSTEM\n")
    g.write("SYMMETRY\n")
    g.write("1\n")
    g.write("CELL\n")
    g.write("%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" % (cell[0]/b, cell[1]/cell[0], cell[2]/cell[0],cos(cell[3]),cos(cell[4]),cos(cell[5])))
    g.write("CUTOFF\n")
    g.write(" 220.0\n")
    g.write("&END\n")
    g.close()

    f=open(file,'w')
    f.write("#BEGIN ID=1 FRAME=1 SERIES=f ENERGY=0\n")
    f.write("#SETPBC=%s\n"%sysfile)
    for i,c in enumerate(coords):
        f.write("%4d %2s%8.4f%8.4f%8.4f  %10.3E %10.3E %10.3E\n" %(i+1,types[i],c[0]/b,c[1]/b,c[2]/b,cos(cell[3]),cos(cell[4]),cos(cell[5])))
    f.write("#END\n")
    f.close()

def writebs(file,Acell,coords,types):
    # xbs format 
    f=open(file,'w')
    for i,l in enumerate(coords):
        f.write("atom %s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))  # xyz format
    f.write('\n')
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[0][0],Acell[0][1],Acell[0][2]))
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[1][0],Acell[1][1],Acell[1][2]))
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[2][0],Acell[2][1],Acell[2][2]))
    f.write('\n')
    f.write('spec O  0.30 red\n')
    f.write('spec C  0.25 black\n')
    f.write('spec B  0.25 0.5\n')
    f.write('spec N  0.2  blue\n')
    f.write('spec P  0.35 green\n')
    f.write('spec H  0.15 yellow\n')
    f.write('\n')
    f.write('bonds C O 0 1.6 0.04 black\n')
    f.write('bonds C H 0 1.6 0.04 black\n')
    f.write('bonds O H 0 1.6 0.04 black\n')
    f.write('bonds O O 0 1.7 0.04 black\n')
    f.write('bonds H H 0 1.3 0.04 black\n')
    f.write('bonds C C 0 1.7 0.04 black\n')
    f.write('bonds B C 0 1.7 0.04 black\n')
    f.write('bonds B B 0 1.7 0.04 black\n')
    f.write('bonds N C 0 1.7 0.04 black\n')
    f.close()

def writeobj(file,coords):
    f=open(file,'w')
    f.write("o crystal\n")
    for c in coords:
      f.write("v %f %f %f\n" % (c[0],c[1],c[2]))
    for i,u in enumerate(coords):
      for j,v in enumerate(coords):
        if i!=j: 
          dx=sqrt(sum((u[k]-v[k])**2 for k in range(3)))
          if dx<1.7: f.write("l %d %d\n" % (i+1,j+1))
    f.close()

if __name__ == "__main__":
    sys.exit(main())
