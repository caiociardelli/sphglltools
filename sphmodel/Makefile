# Makefile for SphModel

# Compiler and linker
CC := gcc
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
	bkmns2dd \
	bkmns2ll \
	bkmns2pf \
	bkmns2mean \
	powspec \
	getinfo

# First list of objects
OBJECTS1 := \
	exmath.o \
	legendre.o \
	coordinates.o \
	spharm.o \
	io.o

# Second list of objects
OBJECTS2 := \
	exmath.o \
	legendre.o \
	coordinates.o \
	spharm.o \
	progress.o \
	io.o

# Complete path for binaries and objects
DFT   := $(patsubst %, $(BIN)/%, $(DEFAULT))
OBJS1 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS1))
OBJS2 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS2))

# Command used for cleaning
RM := rm -rf
 
#
# Compilation and linking
#
all: objDirectory binDirectory $(DFT)
	@ echo 'Finished building binary!'

$(BIN)/bkmns2dd: $(OBJS1) $(OBJ)/bkmns2dd.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/bkmns2ll: $(OBJS1) $(OBJ)/bkmns2ll.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/bkmns2pf: $(OBJS1) $(OBJ)/bkmns2pf.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/bkmns2mean: $(OBJS2) $(OBJ)/bkmns2mean.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/powspec: $(OBJS1) $(OBJ)/powspec.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(BIN)/getinfo: $(OBJ)/getinfo.o
	@ echo 'Building binary using $(GCC) linker: $@'
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIB)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(OBJ)/exmath.o: $(SHARED)/exmath.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/exmath.c -o $@
	@ echo ' '

$(OBJ)/legendre.o: $(SHARED)/legendre.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/legendre.c -o $@
	@ echo ' '

$(OBJ)/coordinates.o: $(SHARED)/coordinates.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/coordinates.c -o $@
	@ echo ' '

$(OBJ)/spharm.o: $(SHARED)/spharm.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/spharm.c -o $@
	@ echo ' '

$(OBJ)/io.o: $(SHARED)/io.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/io.c -o $@
	@ echo ' '

$(OBJ)/progress.o: $(SHARED)/progress.c $(SETUP)
	@ echo 'Building target using $(MPICC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SHARED)/progress.c -o $@
	@ echo ' '

$(OBJ)/bkmns2dd.o: $(SRC)/bkmns2dd.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/bkmns2dd.c -o $@
	@ echo ' '

$(OBJ)/bkmns2ll.o: $(SRC)/bkmns2ll.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/bkmns2ll.c -o $@
	@ echo ' '

$(OBJ)/bkmns2pf.o: $(SRC)/bkmns2pf.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/bkmns2pf.c -o $@
	@ echo ' '

$(OBJ)/bkmns2mean.o: $(SRC)/bkmns2mean.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/bkmns2mean.c -o $@
	@ echo ' '

$(OBJ)/powspec.o: $(SRC)/powspec.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/powspec.c -o $@
	@ echo ' '

$(OBJ)/getinfo.o: $(SRC)/getinfo.c $(SETUP)
	@ echo 'Building target using $(GCC) compiler: $@'
	$(CC) $(FLAGS) $(INC) -c $(SRC)/getinfo.c -o $@
	@ echo ' '

objDirectory:
	@ mkdir -p $(OBJ)

binDirectory:
	@ mkdir -p $(BIN)

clean:
	$(RM) $(OBJ)/ $(BIN)/ *.dat *.cpt *.xyp *.history *.grd *.pdf

.PHONY: all clean
