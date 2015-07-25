#include <stdio.h>
#include <math.h>

#ifndef __mie_h__
#define __mie_h__ 1

const double pi = 3.1415926535897932384;
const double m = 0.8; // Index of refraction (used in Mie scattering)
const double r = 0.00001; // Radius of spherical scattering particles (Mie scattering)
const double lambda = 500; // Wavelength of incident light
const double alpha = (2*pi*r) / lambda; // Size parameter for scattering
const int N = 2; // Number of iterations for accuracy of mie.
const int nangles = 180; // Angular resolution of mie scattering

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

double Legendre_poly(int iter, double x) {
    switch (iter) { 
        case 0:
            return 1;
        case 1:
            return x;
        case 2:
            return 0.5*(3*x*x-1);
        default:
            printf("ERROR: Legendre polynomial called beyond approximation limit.\n");
            return -1;
    }
}

double Legendre_poly_deriv(int iter, double x) {
    switch (iter) { 
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 3*x;
        default:
            printf("ERROR: Legendre polynomial called beyond approximation limit.\n");
            return -1;
    }
}

double a_n(int iter, double x){
    return (bessel1_prime(iter, m*x)*bessel1(iter,x) - m*bessel1(iter,m*x)*bessel1_prime(iter,x)) / 
            (bessel1_prime(iter,m*x)*bessel2(iter,x) - m*bessel1(iter,m*x)*bessel2_prime(iter,x));
}

double b_n(int iter, double x){
    return (m*bessel1_prime(iter, m*x)*bessel1(iter,x) - bessel1(iter,m*x)*bessel1_prime(iter,x)) / 
            (m*bessel1_prime(iter,m*x)*bessel2(iter,x) - bessel1(iter,m*x)*bessel2_prime(iter,x));
}

double mie_pi(int iter, double x) { return (1/sin(x))*Legendre_poly(iter, x); }
double mie_tau(int iter, double x) { return Legendre_poly_deriv(iter, x); }

double scatter_magnitude_1(double theta) {
    double mag = 0;
    for (int i = 0; i < N; i += 1) {
        double term = (a_n(i,alpha)*mie_pi(i,theta)+b_n(i,alpha)*mie_tau(i,theta)); 
        if (i == 0) { mag += term;}
        else { mag += term*((2*i+1)/(i*(i+1))); }
    }
    return mag;
}

double scatter_magnitude_2(double theta) {
    double mag = 0;
    for (int i = 0; i < N; i += 1) {
        double term = (b_n(i,alpha)*mie_pi(i,theta)+a_n(i,alpha)*mie_tau(i,theta)); 
        if (i == 0) { mag += term;}
        else { mag += term*((2*i+1)/(i*(i+1))); }
    }
    return mag;
}

double* mie_data_s11(void) {
    double* result_array = new double[nangles];    
    for (int i = 0; i < nangles; i += 1) {
        double S2 = scatter_magnitude_2(i);
        double S1 = scatter_magnitude_1(i);
        result_array[i] = 0.5*(S2*S2 + S1*S1);
    }
    return result_array;
}

double* mie_data_s12(void) {
    double* result_array = new double[nangles];    
    for (int i = 0; i < nangles; i += 1) {
        double S2 = scatter_magnitude_2(i);
        double S1 = scatter_magnitude_1(i);
        result_array[i] = 0.5*(S2*S2 - S1*S1);
    }
    return result_array;
}

#endif
