SHELL:=/bin/bash

OBJECTS = main.o MIEV0.o ErrPack.o 

all:
	gfortran -c *.f
	g++ -c main.cpp
	#gfortran -o main *.o -lstdc++ -lm
	g++ -O ${OBJECTS} -lm -lgfortran -o main

clean:
	rm *.o
