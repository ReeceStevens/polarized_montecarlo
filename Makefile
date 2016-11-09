SHELL:=/bin/bash

OBJECTS = main.o MIEV0.o ErrPack.o

FLAGS= -DDEBUG -O0

C_COMPILER=g++
# F_COMPILER=gfortran
F_COMPILER=${C_COMPILER}

all:
	${F_COMPILER} -c *.f
	${C_COMPILER} ${FLAGS} -c main.cpp
	#gfortran -o main *.o -lstdc++ -lm
	${C_COMPILER} ${FLAGS} -O ${OBJECTS} -lm -lgfortran -o main

clean:
	rm *.o main

test:
	make all
	./main
