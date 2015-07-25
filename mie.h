#include <stdio.h>
#include <math.h>


double bessel1(int iter, double size) {
    switch (iter) { 
        case 0:
            return sin(size) / size;
        case 1:
            return (sin(size) / (size*size)) - (cos(size) / size);
        case 2:
            return ((3/(size*size*size)) - (1/size))*sin(size) - (3/(size*size*size))*cos(size);
        default:
            printf("ERROR: bessel function called beyond approximation limit.\n");
            return -1;
    }
}

double bessel2(int iter, double size) {
    switch (iter) { 
        case 0:
            return -cos(size) / size;
        case 1:
            return -(cos(size) / (size*size)) - (sin(size) / size);
        case 2:
            return -((3/(size*size*size)) - (1/size))*cos(size) - (3/(size*size*size))*sin(size);
        default:
            printf("ERROR: bessel function called beyond approximation limit.\n");
            return -1;
    }
}

double bessel1_prime(int iter, double size) {
    switch (iter) {
        case 0:
            return ((1/size)*cos(size)) - (1/(size*size))*(sin(size));
        case 1:
            return (2/(size*size))*cos(size) - (2/(size*size*size))*sin(size) + (1/size)*sin(size);
        case 2:
            return ((-9/(size*size*size*size)) + 4/(size*size))*sin(size) + ((9/(size*size*size)) - (1/size))*cos(size);
        default:
            printf("ERROR: bessel function called beyond approximation limit.\n");
            return -1;
    }
}

double bessel2_prime(int iter, double size) {
    switch (iter) {
        case 0:
            return (1/(size*size))*cos(size) + sin(size)/size;
        case 1:
            return ((2/(size*size*size)) - (1/(size*size)))*cos(size) + (2/(size*size))*sin(size);
        case 2:
            return ((-9/(size*size*size*size)) + (4/(size*size)))*cos(size) - ((9/(size*size*size)) - (1/size))*sin(size);
        default:
            printf("ERROR: bessel function called beyond approximation limit.\n");
            return -1;
    }
}


