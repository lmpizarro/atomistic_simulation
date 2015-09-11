MAKEFILE = Makefile
EXE = ising
FC = gfortran #IFORT # /OPT/INTEL/COMPILER70/IA32/BIN/IFC  

# DIrectorio para los modulos
MOD_DIR = ./mod
# Directorio para los objetos
OBJ_DIR = ./obj
# Lista con todos los directorios
BUILD_DIR = $(OBJ_DIR) $(MOD_DIR)

# Warning: the debugger doesn't get along with the optimization options
# So: not use -O3 WITH -g option
# Optimizacion
FLAGS = -O3
# Que buscque modulos en el directorio adecuado  
#FLAGS += -J$(MOD_DIR)  ! Tambien reemplaza al -I . No entiendo.

# Define file extensions
.SUFFIXES:
.SUFFIXES: .f90 .o  

# Lista de archivos con el orden correcto de acuerdo a sus
# dependencias. Ver https://github.com/ZedThree/fort_depend.py
# para hacerlo de forma automatica. O aprender a usar cmake:
# http://brunoabinader.com/?p=29
OBJS_F = ziggurat.o usozig.o isingmods.o $(EXE).o
# Camino completo en donde estaran los objetos
OBJS = $(patsubst %.o, $(OBJ_DIR)/%.o, $(OBJS_F))
# Compruebo que existan los subdirectorios y compila
all:  checkdirs $(EXE)

$(OBJ_DIR)/%.o: %.f90  
	$(FC) $(FLAGS) -J$(MOD_DIR) -c $< -o $@   

$(EXE): $(OBJS) Makefile 
	$(FC) $(FLAGS) -I$(MOD_DIR) -o $@ $(OBJS) 

$(EXE).o: $(EXE).f90 ziggurat.o usozig.o isingmods.o 

checkdirs : $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir $@ 

install: $(EXE) ./bin
	@ln $(EXE) ./bin
	@cp input.dat histograma.plg graficar.sh ./bin

./bin:
	@mkdir ./bin

clean:
	rm $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod	

distclean:
	rm -fr $(BUILD_DIR)
	rm ./bin/$(EXE)

.PHONY: clean distclean
