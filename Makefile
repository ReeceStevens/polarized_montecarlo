SHELL:=/bin/bash

OBJECTS = main.o MIEV0.o ErrPack.o 

C_COMPILER=g++-5
F_COMPILER=gfortran

all:
	${F_COMPILER} -c *.f
	${C_COMPILER} -c main.cpp
	#gfortran -o main *.o -lstdc++ -lm
	${C_COMPILER} -O ${OBJECTS} -lm -lgfortran -o main

clean:
	rm *.o
	rm main

test:
	make all
	./main
