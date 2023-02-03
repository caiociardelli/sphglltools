# Makefile for SphGLLTools

# Compiler and linker
CC := gcc
# MPI compiler and linker
MPICC := mpicc
# Flags for compiler
FLAGS := \
	-Wall \
	-Wextra \
	-std=c99 \
	-O3
# Source directory
SRC := src
# Shared objects directory
SHARED := $(SRC)/shared
# Configuration files directory
SETUP := setup
# Configuration files
CONFIG := $(wildcard $(SETUP)/*.h)
# Includes
INC := -I$(SETUP) -I$(SRC)/headers
# Objects directory
OBJ := obj
# Libraries
LIB := -lm
# Binaries directory
BIN := bin

# Default targets
DEFAULT := \
	gll2bk \
	gll2dd \
	gll2gll \
	gll2ll \
	gll2mean \
	gll2pf \
	gll2mns \
	setspl \
	smooth \
	getinfo

# First list of objects
OBJECTS1 := \
	exmath.o \
	legendre.o \
	coordinates.o \
	io.o \
	boundaries.o \
	progress.o

# Second list of objects
OBJECTS2 := \
	exmath.o \
	legendre.o \
	coordinates.o \
	io.o \
	stretch.o \
	boundaries.o \
	progress.o

# Third list of objects
OBJECTS3 := \
	exmath.o \
	legendre.o \
	coordinates.o \
	io.o \
	spharm.o \
	boundaries.o \
	progress.o

# Fourth list of objects
OBJECTS4 := \
	exmath.o \
	coordinates.o \
	io_smooth.o \
	boundaries_smooth.o \
	profile.o \
	progress.o

# Complete path for binaries and objects
DFT   := $(patsubst %, $(BIN)/%, $(DEFAULT))
OBJS1 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS1))
OBJS2 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS2))
OBJS3 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS3))
OBJS4 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS4))

# Command used for cleaning
RM := rm -rf

#
# Compilation and linking
#
all: objDirectory binDirectory $(DFT)
	@ echo 'Finished building binary!'

$(BIN)/gll2bk: $(OBJS1) $(OBJ)/gll2bk.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2dd: $(OBJS1) $(OBJ)/gll2dd.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2gll:$(OBJS2) $(OBJ)/gll2gll.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2ll: $(OBJS1) $(OBJ)/gll2ll.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2mean: $(OBJS1) $(OBJ)/gll2mean.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2pf: $(OBJS1) $(OBJ)/gll2pf.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/gll2mns: $(OBJS3) $(OBJ)/gll2mns.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/setspl: $(OBJS3) $(OBJ)/setspl.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/smooth: $(OBJS4) $(OBJ)/smooth.o
	@ echo 'Building binary using $(MPICC) linker: $@'
	$(MPICC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/getinfo: $(OBJ)/getinfo.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(OBJ)/exmath.o: $(SHARED)/exmath.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/exmath.c -o $@
	@ echo ' '

$(OBJ)/legendre.o: $(SHARED)/legendre.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/legendre.c -o $@
	@ echo ' '

$(OBJ)/coordinates.o: $(SHARED)/coordinates.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/coordinates.c -o $@
	@ echo ' '

$(OBJ)/io.o: $(SHARED)/io.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/io.c -o $@
	@ echo ' '

$(OBJ)/io_smooth.o: $(SHARED)/io_smooth.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/io_smooth.c -o $@
	@ echo ' '

$(OBJ)/spharm.o: $(SHARED)/spharm.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/spharm.c -o $@
	@ echo ' '

$(OBJ)/stretch.o: $(SHARED)/stretch.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/stretch.c -o $@
	@ echo ' '

$(OBJ)/boundaries.o: $(SHARED)/boundaries.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/boundaries.c -o $@
	@ echo ' '

$(OBJ)/boundaries_smooth.o: $(SHARED)/boundaries_smooth.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/boundaries_smooth.c -o $@
	@ echo ' '

$(OBJ)/profile.o: $(SHARED)/profile.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/profile.c -o $@
	@ echo ' '

$(OBJ)/progress.o: $(SHARED)/progress.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SHARED)/progress.c -o $@
	@ echo ' '

$(OBJ)/gll2bk.o: $(SRC)/gll2bk.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2bk.c -o $@
	@ echo ' '

$(OBJ)/gll2dd.o: $(SRC)/gll2dd.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2dd.c -o $@
	@ echo ' '

$(OBJ)/gll2gll.o: $(SRC)/gll2gll.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2gll.c -o $@
	@ echo ' '

$(OBJ)/gll2ll.o: $(SRC)/gll2ll.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2ll.c -o $@
	@ echo ' '

$(OBJ)/gll2mean.o: $(SRC)/gll2mean.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2mean.c -o $@
	@ echo ' '

$(OBJ)/gll2pf.o: $(SRC)/gll2pf.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2pf.c -o $@
	@ echo ' '

$(OBJ)/gll2mns.o: $(SRC)/gll2mns.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/gll2mns.c -o $@
	@ echo ' '

$(OBJ)/setspl.o: $(SRC)/setspl.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/setspl.c -o $@
	@ echo ' '

$(OBJ)/smooth.o: $(SRC)/smooth.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(MPICC) $(FLAGS) $(INC) -c $(SRC)/smooth.c -o $@
	@ echo ' '

$(OBJ)/getinfo.o: $(SRC)/getinfo.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/getinfo.c -o $@
	@ echo ' '

objDirectory:
	@ mkdir -p $(OBJ)

binDirectory:
	@ mkdir -p $(BIN)

clean:
	$(RM) $(OBJ)/ $(BIN)/ *.dat *.cpt *.xyp *.history *.grd *.pdf

.PHONY: all clean
