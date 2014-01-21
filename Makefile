# Use Intel C compiler 
CC = icc 
CPP = icpc
# Compile-time flags 
CFLAGS = -std=c99 -lstdc -O3 -fopenmp
CPPFLAGS = -fopenmp
DEBUG = 

all: jdh_solv

jdh_solv: slv.o
	$(CPP) ./src/slv.o $(CPPFLAGS) $(DEBUG) -o jdh_solv

slv.o: 
	$(CC) -c ./src/slv.c $(CFLAGS) $(DEBUG) -o ./src/slv.o
	
clean:
	rm -rf ./src/*o jdh_solv