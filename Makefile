# Use Intel C compiler 
CC = icc 
CPP = icpc
# Compile-time flags 
CFLAGS = -std=c99 -llapack -lblas -lm -lstdc -O3 -openmp -xSSE4.2 -fno-alias
CPPFLAGS =  -llapack -lblas -lm -openmp -O3 -xSSE4.2 -fno-alias
DEBUG = 
#DEBUG = -DDEBUG -DDBGMGM -DDBGMGM

all: jdh_solv

jdh_solv: slv.o mgmres.o
	$(CPP) ./src/slv.o ./src/mgmres.o $(CPPFLAGS) $(DEBUG) -o jdh_solv

slv.o: 
	$(CC) -c ./src/slv.c $(CFLAGS) $(DEBUG) -o ./src/slv.o
	
mgmres.o:
	$(CC) -c ./src/mgmres.c $(CFLAGS) $(DEBUG) -o ./src/mgmres.o
	
clean:
	rm -rf ./src/*o jdh_solv ./*.dat
	
