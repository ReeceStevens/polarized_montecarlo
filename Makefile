SHELL:=/bin/bash

all:
	gfortran -c *.f
	gcc -c main.cpp
	#gfortran -o main *.o -lstdc++ -lm
	gcc -o main main.o MIEV0.o ErrPack.o -lstdc++ -lgfortran -lm 

clean:
	rm *.o
