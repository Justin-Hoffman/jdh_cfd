# Use Intel C compiler 
CC = icc 
CPP = icpc
# Compile-time flags 
#CFLAGS = -std=c99 -llapack -lblas -lm -lstdc -openmp -fno-alias -g -O0 -traceback
CFLAGS = -std=c99 -llapack -lblas -lm -lstdc -openmp -fno-alias -O3
#CPPFLAGS =  -llapack -lblas -lm -openmp -fno-alias -g -O0 -traceback
CPPFLAGS =  -llapack -lblas -lm -openmp -fno-alias -O3
DEBUG = 
#DEBUG = -DDEBUG -DDBGVBE -DDBGMGM -DDBGMEM -DDBGBCS -DDBGLVLST

all: jdh_solv

jdh_solv: slv.o lvlset.o mgmres.o util.o
	$(CPP) ./src/slv.o ./src/lvlset.o ./src/mgmres.o ./src/util.o $(CPPFLAGS) $(DEBUG) -o jdh_solv

slv.o: 
	$(CC) -c ./src/slv.c $(CFLAGS) $(DEBUG) -o ./src/slv.o
	
lvlset.o: 
	$(CC) -c ./src/lvlset.c $(CFLAGS) $(DEBUG) -o ./src/lvlset.o
	
mgmres.o:
	$(CC) -c ./src/mgmres.c $(CFLAGS) $(DEBUG) -o ./src/mgmres.o
	
util.o:
	$(CC) -c ./src/util.c $(CFLAGS) $(DEBUG) -o ./src/util.o
	
	
clean:
	rm -rf ./src/*o jdh_solv ./*.dat
	
