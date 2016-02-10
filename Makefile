all:
	rm *.o
	gfortran -c *.f
	gcc -c main.cpp
	#gfortran -o main *.o -lstdc++ -lm
	gcc -o main ErrPack.o MIEV0.o main.o -lstdc++ -lgfortran -lm 
