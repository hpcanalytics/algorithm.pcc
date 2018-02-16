MPICXX=g++
CPPFLAGS=-DX86 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES=1
#CPPFLAGS=-Wall
LDFLAGS=-lmpi -lmpi++

SRC = ./src
BIN = ../../Bin

parcc:utility.hpp Set.hpp ParClusterCo.hpp ParClusterCo.cpp
	module load mvapich2/gcc;	mpicxx -o parcc ParClusterCo.cpp
	
clean:
	rm -f $(BIN)/parcc	
