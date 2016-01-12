#!/usr/local/bin/python
#Authors: Laalitha Liyanage, Sungho Kim
#Purpose: Calculate material properties of fcc system
# https://icme.hpc.msstate.edu/mediawiki/index.php/M.E.A.M_Python_Script


import os
import re
import sys
import math
import numpy as np
#========================================================================================================================================
#--------------------------------------------------Default setting-----------------------------------------------------------------------
#========================================================================================================================================

datafile="datafile"

# Default lattice parameter
a=4.05 

#Paramters of interatomic potential
pair_style = "meam"
pair_coeff = "* * meamf AlS meam.alsimgcufe AlS"

#Constants to convert units
unit_conversion1 = 160.217646 # 1 eV/A^3 = 1.60217646E-19 / 1E-30 Pa = 160.217646 GPa
unit_conversion2 = 16021.765 #1 eV/A^3 to mJ/m^2

#Command to execute LAMMPS
executable = "mpirun -np 8 ~/packages/lammps-10Aug15/src/lmp_mpi"
#executable = "qsub run.qsub"

#Percentage strain and no. of points for elastic constant calculation
perc = .1
nt = 6

#Length of vacuum for surface energy calculation
vacuum =  15.0

#=============================================================================================================================#
def get_fcc_coords(a,nx=1,ny=1,nz=1,vac=0):#Generate coordinates for of FCC
 
  xa = []; ya = []; za = []; ty = []
  ax = a
  ay = a
  az = a
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):        
        if vac == 0:
                xa.append( (0.0+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)
                xa.append( (0.0+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                xa.append( (0.5+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                xa.append( (0.5+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)
        if vac == 1:
                if i==nx/2 and j==ny/2 and k==nz/2:
                        xa.append( (0.0+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                        xa.append( (0.5+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                        xa.append( (0.5+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)

                else:
                        xa.append( (0.0+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)
                        xa.append( (0.0+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                        xa.append( (0.5+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
                        xa.append( (0.5+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)

  return ty,xa,ya,za,ax*nx,ay*ny,az*nz

def distort_simbox(xa,ya,za,box,strain_tensor=[0,0,0,0,0,0]):#Apply strain and calculate the 9 triclinic box parameters for LAMMPS

  #Initial orthorgonal box edges taken as lattice translation vectors
  ax,ay,az = box[0],box[1],box[2]

  #Strains
  e1,e2,e3,e4,e5,e6 = strain_tensor

  #Application of strain to the lattice vectors
  ax_p = [(1.0+e1)*ax,0.5*e6*ax,0.5*e5*ax]
  ay_p = [0.5*e6*ay,(1.0+e2)*ay,0.5*e4*ay]
  az_p = [0.5*e5*az,0.5*e4*az,(1.0+e3)*az]
   
  #The new lattice vectors (cell vectors)
  cell =  [ax_p,ay_p,az_p]
  
  #Get length of cell vectors
  a = np.linalg.norm(cell[0]) 
  b = np.linalg.norm(cell[1]) 
  c = np.linalg.norm(cell[2]) 
  
  #Calculate angles between cell vectors
  alpha = np.arccos( np.vdot(cell[1], cell[2]) / (b * c) )
  beta = np.arccos( np.vdot(cell[0], cell[2]) / (a * c) )
  gamma = np.arccos( np.vdot(cell[0], cell[1]) / (a * b) )
 
  # a_LAMMPS = (xhi-xlo,0,0); b_LAMMPS = (xy,yhi-ylo,0); c_LAMMPS = (xz,yz,zhi-zlo)
  xlo = ylo = zlo = 0.0
  # this choice of origin simplifies things:
  # a_LAMMPS = (xhi,0,0); b_LAMMPS = (xy,yhi,0); c_LAMMPS = (xz,yz,zhi)
  
  xhi = a # a_LAMMPS
  xy = np.cos(gamma) * b# b_LAMMPS
  yhi = np.sin(gamma) * b
  xz = np.cos(beta) * c# c_LAMMPS
  yz = ( b * c * np.cos(alpha) - xy * xz ) / yhi
  zhi = np.sqrt( c**2 - xz**2 - yz**2 )
 
  #LAMMPS vectors
  
  cell_lammps = np.array([[xhi-xlo,0,0],[xy,yhi-ylo,0],[xz,yz,zhi-zlo]])

  bx = [xhi-xlo,0,0]
  by = [xy,yhi-ylo,0]
  bz = [xz,yz,zhi-zlo]
 
  #print cell, cell_lammps
  
  # IMPORTANT: need vector-(rotation-)matrix product (instead of matrix-vector product) here,
  #            cell vectors are ROW VECTORS (also see above)
  
  rotation = np.dot(np.linalg.inv(cell), cell_lammps)

  #print rotation

  for i in range(len(xa)):
    r = [xa[i],ya[i],za[i]]
    [x,y,z] = np.dot(r,rotation)
    #print r, [x,y,z]
    xa[i]=x;ya[i]=y;za[i]=z

  return xa,ya,za,bx,by,bz

def gen_datafile(ty,xa,ya,za,bx,by,bz):
  fout = open(datafile,"w")
  fout.write("cementite a = %f, b = %f, c = %f\n\n" % (bx[0],by[1],bz[2]))
  fout.write("%d atoms\n"%len(xa))
  fout.write("1 atom types\n")
  fout.write(" 0.0  %22.16f   xlo xhi\n"%bx[0])
  fout.write(" 0.0  %22.16f   ylo yhi\n"%by[1])
  fout.write(" 0.0  %22.16f   zlo zhi\n"%bz[2])
  fout.write(" %22.16f  %22.16f %22.16f xy xz yz\n"%(by[0],bz[0],bz[1]))
  fout.write("\nAtoms\n\n")
  for i in range(len(xa)):
    fout.write("%4d %3d %22.16f %22.16f %22.16f\n"%(i+1,ty[i],xa[i],ya[i],za[i]))
  fout.close()
  return len(xa)

def gen_infile(relax=0,box_relax=0):# Create input command script for lammps
  if relax == 1:
      fname = "infile_relax"
  else:
      fname = "infile_static"

  fout = open(fname,'w')
  fout.write("units           metal\n")
  fout.write("boundary        p p p\n")
  fout.write("atom_style      atomic\n")
  fout.write("read_data       datafile\n")
  fout.write("pair_style      %s\n"%pair_style)
  fout.write("pair_coeff      %s\n"%pair_coeff)
  fout.write("neighbor        2.0 bin\n")
  fout.write("neigh_modify    delay 10\n")
  fout.write("dump            1 all custom 1 dump id type x y z\n")
  fout.write("log             log\n")
  fout.write("thermo_style    custom step atoms pe ke etotal press vol lx ly lz\n")
  fout.write("thermo_modify   line one format float %.16g\n")
  
  if relax ==0 and box_relax ==0:
      fout.write("run             0\n")
  elif relax == 1 and box_relax == 0:
      fout.write("min_style       cg\n")
      fout.write("min_modify      line quadratic\n")
      fout.write("minimize        1.0e-30 1.0e-20 100000 1000000\n")
  elif relax == 1 and box_relax ==1:
      fout.write("fix             1 all box/relax iso 0.0 vmax 0.001\n")
      fout.write("min_style       cg\n")
      fout.write("min_modify      line quadratic\n")
      fout.write("minimize        1.0e-30 1.0e-20 100000 1000000\n")
  fout.close()
  return fname

def get_field_from_log_lammps(field):# Get data from Lammps log file
  os.system("grep Loop log -B1|grep -v Loop|awk '{print $%s}' > eout"%(field))
  fin = file("eout",'r')
  line = fin.readline().split()
  fin.close()
  os.system("rm eout")
  return (float(line[0]))

def get_E0():#Relax structure and get equilibrium energy and volume.
  ty,xa,ya,za,ax,ay,az = get_fcc_coords(a)
  xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[0,0,0,0,0,0])
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=1)))

  N =  get_field_from_log_lammps("2")
  E0 = get_field_from_log_lammps("3")
  V0 = get_field_from_log_lammps("7")
  a0 = get_field_from_log_lammps("8")
  os.system("mv log log.E0")
  os.system("mv dump dump.E0")
  return(N,E0, V0, a0)

##---------------Function definitions relating to Elastic Constant calculation-----------------------------------##

def get_cxx(cname,a0,a,b,c):
  energy = []
  os.system("echo %s > report.lammps"%(cname))
  t0,t1 = -perc/100,perc/100
  dt = (t1 - t0)/nt

  ty,xa,ya,za,ax,ay,az = get_fcc_coords(a0)

  for i in range(nt+1):

    t = t0+i*dt

    #Generate lattice according to strain tensor
    if cname == 'c11': 
      xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[t,0,0,0,0,0])
    elif cname == 'c12':
      xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[t,-t,0,0,0,0])
    elif cname == 'c44':
      xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[0,0,0,2*t,0,0])

    gen_datafile(ty,xa,ya,za,bx,by,bz)
    #Run LAMMPS and get energy
    os.system("echo 't=%f' >> report.lammps "%(t))
    os.system("%s -in %s >> report.lammps"%(executable,gen_infile(relax=1,box_relax=0)))
    energy.append(get_field_from_log_lammps("5"))

  fout = open("summary","w")
  for i in range(nt+1):
    fout.write("%22.16f %22.16f\n"%(t0+i*dt,energy[i]))
  fout.close()
  gpfile = gen_gpfile_for_parabola_fit(a*a0**3,b,c,cname)
  os.system("gnuplot %s > gnuplot.report 2>&1"%gpfile)
  a,b,c = get_fitted_param(a*a0**3,b,c,"fit.log")
  os.system("mv fit.log %s_fit.log"%cname)
  os.system("mv report.lammps %s_report.lammps"%cname)
  return a/a0**3

def gen_gpfile_for_parabola_fit(a,b,c,prefix):
#Generate input script for GNUplot for fitting
  os.system("rm -f fit.log")
  fn = "%s_fit.gp"%prefix
  fout = open(fn,"w")
  fout.write('set term gif\n')
  fout.write('set output "%s_fit.gif"\n'%prefix)
  fout.write("f(x) = a*x**2 + b*x+ c\n")
  b = 0.0
  fout.write("a = %f; b = %f; c = %f\n"%(a,b,c))
  fout.write("FIT_LIMIT = 1e-9\n")
  if b == 0.0:
    fout.write("fit f(x) 'summary' via a,c\n")
  else:
    fout.write("fit f(x) 'summary' via a,b,c\n")
  fout.write('fx_title = sprintf("(%f)*x**2+(%f)*x+(%f)",a,b,c)\n')
  fout.write("plot 'summary' t 'MEAM' with linespoints 1 5,f(x) t fx_title")
  #fout.write("plot 'summary' with linespoints 1 5,f(x)")
  fout.close()
  return fn

def get_fitted_param(a,b,c,fn):
#Extract fitted parameters from GNUplot output
  fin = file(fn,"r")
  while True:
    line = fin.readline()
    if line == "": sys.exit(123)
    if re.search("^Final set of parameters",line): break
  for i in range(2): line = fin.readline()
  for i in range(2):
    line = fin.readline().split()
    if   len(line)>0 and line[0] == "a": a = float(line[2])
    elif len(line)>0 and line[0] == "b": b = float(line[2])
    elif len(line)>0 and line[0] == "c": c = float(line[2])
  fin.close()
  return a,b,c

##---------------Function definitions relating to Surface Energy calculation-------------------------------##

#--------------------------Surface (100)---------------------------------------------------------

def gen_coords_fcc_100(a,nx=1,ny=1,nz=1):#Generate coordinates for fcc (100) surface structure
  xa=[]; ya=[]; za=[];ty = []
  bx = [0]*3;by = [0]*3; bz = [0]*3

  bx[0],by[1],bz[2] = a*nx,a*ny,a*nz+vacuum
  x,y,z = bx,by,bz

  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(0 + i*a); ya.append(  0 + j*a); za.append(  0 + k*a); ty.append(1)
        xa.append(  0 + i*a); ya.append(a/2 + j*a); za.append(a/2 + k*a);ty.append(1)
        xa.append(a/2 + i*a); ya.append(  0 + j*a); za.append(a/2 + k*a);ty.append(1)
        xa.append(a/2 + i*a); ya.append(a/2 + j*a); za.append(  0 + k*a);ty.append(1)
  return ty,xa,ya,za,bx,by,bz

def get_surf_energy_100(a0,E0):#Calculate fcc (100) surface formation energy
  ty,xa,ya,za,bx,by,bz = gen_coords_fcc_100(a0,2,2,4)
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=0)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e100 = (E-N*E0)/2/(get_field_from_log_lammps("8")*get_field_from_log_lammps("9"))
  os.system("mv log log.surf100")
  os.system("mv dump dump.surf100")
  return(e100)

#--------------------------Surface (110)--------------------------------------------------------------

def gen_coords_fcc_110(a,nx=1,ny=1,nz=1):#Generate coordinates for fcc (110) surface structure
  xa=[]; ya=[]; za=[]; ty=[]
  bx = [0]*3;by = [0]*3; bz = [0]*3

  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4. * a
  y2 = math.sqrt(6)/4. * a
  y3 = math.sqrt(6)/6. * a
  y4 = math.sqrt(6)*5./12. * a
  y5 = math.sqrt(6)*2./6. * a
  y6 = math.sqrt(6)/12 * a
  z3 = math.sqrt(3)/3. * a
  z5 = math.sqrt(3)*2./3. * a
  bx[0],by[1],bz[2] = ax*nx + vacuum, ay*ny, az*nz
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(x0+k*az);ty.append(1)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(x0+k*az);ty.append(1)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(z3+k*az);ty.append(1)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(z3+k*az);ty.append(1)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(z5+k*az);ty.append(1)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(z5+k*az);ty.append(1)
  
  return ty,xa,ya,za,bx,by,bz

def get_surf_energy_110(a0,E0):#Calculate fcc (110) surface formation energy
  ty,xa,ya,za,bx,by,bz = gen_coords_fcc_110(a0,4,2,2)
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=0)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e110 = (E-N*E0)/2/(get_field_from_log_lammps("9")*get_field_from_log_lammps("10"))
  os.system("mv log log.surf110")
  os.system("mv dump dump.surf110")
  return e110

#--------------------------Surface (111)-----------------------------------------------------------------------

def gen_coords_fcc_111(a,nx=1,ny=1,nz=1):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  xa=[]; ya=[]; za=[];ty =[]
  bx = [0]*3;by = [0]*3; bz = [0]*3

  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4 * a
  y2 = math.sqrt(6)/4 * a
  y3 = math.sqrt(6)/6 * a
  y4 = math.sqrt(6)*5/12 * a
  y5 = math.sqrt(6)*2/6 * a
  y6 = math.sqrt(6)/12 * a
  bx[0],by[1],bz[2] = ax*nx, ay*ny, az*nz+vacuum
  for i in range(nx):
    for j in range(ny):
      layer = 0
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)

  return ty,xa,ya,za,bx,by,bz

def get_surf_energy_111(a0,E0):#Calculate fcc (111) surface formation energy
  ty,xa,ya,za,bx,by,bz = gen_coords_fcc_111(a0,2,2,4)
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=0)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e111 = (E-N*E0)/2/(get_field_from_log_lammps("8")*get_field_from_log_lammps("9"))
  os.system("mv log log.surf111")
  os.system("mv dump dump.surf111")
  return e111

##---------------Function definitions relating to Generalized Stacking Fault Energy------------------------------------##

def gen_data_for_stacking_fault(a,nx=1,ny=1,type="I"):#Generate coordinates for FCC surface: 110:x, 112:y, 111:z
  
  xa=[]; ya=[]; za=[];ty = []
  bx = [0]*3;by = [0]*3; bz = [0]*3

  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4 * a
  y2 = math.sqrt(6)/4 * a
  y3 = math.sqrt(6)/6 * a
  y4 = math.sqrt(6)*5/12 * a
  y5 = math.sqrt(6)*2/6 * a
  y6 = math.sqrt(6)/12 * a
  if type == "I":
    nlayer = 11
    nfaults = 1
  elif type == "E":
    nlayer = 13
    nfaults = 1

  bx[0],by[1],bz[2] = ax*nx, ay*ny, az/3.0*nlayer

  for i in range(nx):
    for j in range(ny):
      layer = 0
      for k in range(2):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
      if type == "I":
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1) 
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
      elif type == "E":                                                           
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
        for k in range(2):
          xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
          xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
          xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
          xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
          xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)            ;ty.append(1)
          xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1;ty.append(1)
  return ty,xa,ya,za,bx,by,bz,nfaults

def get_stacking_fault(a0,E0,type='I'):
  ty,xa,ya,za,bx,by,bz,nfaults = gen_data_for_stacking_fault(a0,2,2,type)
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=0)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  gsfe = (E-N*E0)/(get_field_from_log_lammps("8")*get_field_from_log_lammps("9"))
  os.system("mv log log.gsfe.%s"%type)
  os.system("mv dump dump.gsfe.%s"%type)
  return gsfe

##-------------------Calculation of vacancy formation energy-----------------------------------------------------##

def get_vac_form_energy(a0,E0):
  ty,xa,ya,za,ax,ay,az = get_fcc_coords(a,4,4,4,vac=1)
  xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[0,0,0,0,0,0])
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=1)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e_vac = E-E0*N
  os.system("mv log log.vac")
  os.system("mv dump dump.vac")
  return(e_vac)
  
##-------------------------Calculation of interstitial energy-----------------------------------------------------##

def get_fcc_coords_with_int(a,nx=1,ny=1,nz=1,int='oct'):#Generate coordinates for of FCC with interstitial atom
 
  xa = []; ya = []; za = []; ty = []
  ax = a
  ay = a
  az = a
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        if int == 'oct' and i==nx/2 and j==ny/2 and k==nz/2:
           xa.append( (0.5+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)  
        if int == 'tet' and i==nx/2 and j==ny/2 and k==nz/2:
           xa.append( (0.25+i)*ax ); ya.append( (0.25+j)*ay ); za.append( (0.25+k)*az ); ty.append(1)  
        xa.append( (0.0+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)
        xa.append( (0.0+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
        xa.append( (0.5+i)*ax ); ya.append( (0.0+j)*ay ); za.append( (0.5+k)*az ); ty.append(1)
        xa.append( (0.5+i)*ax ); ya.append( (0.5+j)*ay ); za.append( (0.0+k)*az ); ty.append(1)
  return ty,xa,ya,za,ax*nx,ay*ny,az*nz

def get_int_energy(a0,E0):
#Octahedral interstitial energy calculation
  ty,xa,ya,za,ax,ay,az = get_fcc_coords_with_int(a,4,4,4,int='oct')
  xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[0,0,0,0,0,0])
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=1)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e_oct = E-E0*N
  os.system("mv log log.oct")
  os.system("mv dump dump.oct")
#Tetrahedral interstitial energy calculation
  ty,xa,ya,za,ax,ay,az = get_fcc_coords_with_int(a,4,4,4,int='tet')
  xa,ya,za,bx,by,bz = distort_simbox(xa,ya,za,[ax,ay,az],strain_tensor=[0,0,0,0,0,0])
  gen_datafile(ty,xa,ya,za,bx,by,bz)
  os.system("%s -in %s > report.lammps"%(executable,gen_infile(relax=1,box_relax=1)))
  N =  get_field_from_log_lammps("2")
  E = get_field_from_log_lammps("3")
  e_tet = E-E0*N
  os.system("mv log log.tet")
  os.system("mv dump dump.tet")
  return(e_oct, e_tet)

#Main program

#Get equilibrium crystal structure data
N,E0,V,a0 = get_E0()
print N, E0, V, a0
#Get elastic constants
c11 = get_cxx("c11",a0,0.5,1.0,E0)*2*unit_conversion1
k12 = get_cxx("c12",a0,0.5,1.0,E0)*2*unit_conversion1
c22=c11
c12 = (c11+c22-k12)/2
c44 = get_cxx("c44",a0,2.0,1.0,E0)/2*unit_conversion1
B=(c11+2*c12)/3
#Get surface energies
E100 = get_surf_energy_100(a0,E0/N)*unit_conversion2
E110 = get_surf_energy_110(a0,E0/N)*unit_conversion2
E111 = get_surf_energy_111(a0,E0/N)*unit_conversion2
#Get stacking fault energies
gsfe_i = get_stacking_fault(a0,E0/N,'I')*unit_conversion2
gsfe_e = get_stacking_fault(a0,E0/N,'E')*unit_conversion2
#Get defect energies
E_vac =  get_vac_form_energy(a0,E0/N)
E_oct,E_tet = get_int_energy(a0,E0/N)

#Write results to file

fout = open("responses.txt",'w')

S = "E_coh\t\t%f\t(-3.35)\teV/atom\n" % (E0/N);fout.write(S)
S = "a_0\t\t%f\t(4.05)\tAngstrom(1e-10m)\n" % (a0);fout.write(S)
S = "V_0\t\t%f\t(16.61)\tAngstrom^3/atom\n" % (V/N);fout.write(S)
S = "C11\t\t%f\t(107)\tGPa\n" % (c11);fout.write(S)
S = "C12\t\t%f\t(61)\tGPa\n" % (c12);fout.write(S)
S = "C44\t\t%f\t(28)\tGPa\n" % (c44);fout.write(S)
S = "B\t\t%f\t(82)\tGPa\n" % (B);fout.write(S)
S = "E_100\t\t%f\t(890)\tmJ/m^2\n" % (E100);fout.write(S)
S = "E_110\t\t%f\t(960)\tmJ/m^2\n" % (E110);fout.write(S)
S = "E_111\t\t%f\t(780)\tmJ/m^2\n" % (E111);fout.write(S)
S = "GSFE_I\t\t%f\t(133)\tmJ/m^2\n" % (gsfe_i);fout.write(S)
S = "GSFE_I\t\t%f\t(133)\tmJ/m^2\n" % (gsfe_e);fout.write(S)
S = "E_vac\t\t%f\t(0.5)\teV\n" % (E_vac);fout.write(S)
S = "E_Oct\t\t%f\t(2.8)\teV\n" % (E_oct);fout.write(S)
S = "E_Tetra\t\t%f\t(3.3)\teV\n" % (E_tet);fout.write(S)

fout.close()
