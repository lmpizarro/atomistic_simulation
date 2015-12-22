# este Makefile queda como ejemplo de referencia
MAKEFILE = Makefile
exe = main 
fcomp = gfortran #ifort # /opt/intel/compiler70/ia32/bin/ifc  
# Warning: the debugger doesn't get along with the optimization options
# So: not use -O3 WITH -g option
flags =  -O3  
# Remote compilation
OBJS = main.o read_param.o combinacion.o globales.o types.o

.SUFFIXES:            # this deletes the default suffixes 
.SUFFIXES: .f90 .o    # this defines the extensions I want 

.f90.o:  
	$(fcomp) -c $(flags) $< 
        

$(exe):  $(OBJS) Makefile 
	$(fcomp) $(flags) -o $(exe) $(OBJS) 


clean:
	rm ./*.o ./*.mod main	

main.o: read_param.o combinacion.o globales.o types.o
read_param.o: read_param.f90 globales.o types.o
types.o: types.f90
globales.o: globales.f90 types.o
combinacion.o: combinacion.f90  globales.o types.o
