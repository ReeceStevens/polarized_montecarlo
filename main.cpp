/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double nphotons = 1000000;
/* Create Global Variables */
// Stokes vector collector for reflected photons
double I_R;	
double Q_R;	
double U_R;	
double V_R;	

// Stokes vector collector for transmitted photons
double I_T;	
double Q_T;	
double U_T;	
double V_T;	

/* Mie Setup */

struct fortran_complex{ float r,i; };
double complx_abs(fortran_complex a) {
    double abs = sqrt(a.r*a.r + a.i*a.i);
    return abs;  
}
fortran_complex complx_mult(fortran_complex a, fortran_complex b) {
    fortran_complex product;
    product.r = a.r*b.r - a.i*b.i;
    product.i = a.r*b.i + a.i*b.r;
    return product; 
}
fortran_complex complx_conj(fortran_complex a) {
    fortran_complex result;
    result.r = a.r;
    result.i = -a.i;
    return result;
}

float pi = 3.1415926536;
double radius = 1.00; // TODO: check units of radius. I think microns?
double wavelength = 0.600; // TODO: microns?  
int nangles = 1000;

double* s11 = new double[nangles];
double* s12 = new double[nangles];
double* s33 = new double[nangles];
double* s43 = new double[nangles];

extern "C" void MIEV0(float* XX, fortran_complex* CREFIN, int* PERFCT, float* MIMCUT, int* ANYANG, 
                    int* NUMANG, 
                   float* XMU, int* NMOM, int* IPOLZN, int* MOMDIM, int* PRNT, float* QEXT, 
			       float* QSCA, float* GQSC, float* PMOM, fortran_complex* SFORW, fortran_complex* SBACK,
                   fortran_complex* S1, fortran_complex* S2, 
                   fortran_complex* TFORW, fortran_complex* TBACK, float* SPIKE );

#include "photons.h"

int main() {

srand(time(NULL));
// Mie scattering variables. Passed into Wiscombe's mie scattering subroutine.
// Initial values taken from Ramella et. al. model.
// INPUT
float XX= (float) (2 * pi * radius / wavelength); // Size parameter
//float XX= 10.0; // Test value
//Complex CREFIN((1.59/1.33),0.0); // fortran_complex refractive index
fortran_complex CREFIN;
CREFIN.r = (1.59/1.33);
CREFIN.i = 0.0;
int PERFCT = 0; // false, refractive index is not infinite.
float MIMCUT = 0.0; // TODO: is this correct for cutoff of imaginary index, since it's already zero?
int ANYANG = 0; // Angles are monotone increasing and mirror symmetric around 90 degrees
int NUMANG = nangles; // number of angles to be evaluated for S1 and S2.
float* XMU = new float[nangles];
// fill out angle array
for(int i = 0; i <nangles; i ++) {
	XMU[i] = (float) cos(pi * i / nangles);
}
int NMOM = 0; // Avoid calculation of PMOM
int IPOLZN = -1; // Not used.
int MOMDIM = 1; // Not used.
int* PRNT = new int[2];
PRNT[0] = 0;
PRNT[1] = 0;
// OUTPUT
float QEXT, QSCA, GQSC; // Extinction efficiency factor, scattering efficiency factor, and asymmetry factor
// Mie scattering amplitudes (what we need!)
fortran_complex* S1 = new fortran_complex[nangles];
fortran_complex* S2 = new fortran_complex[nangles];
//fortran_complex SFORW(0,0); // Forward scattering amplitude S1 at 0 degrees.
//fortran_complex SBACK(0,0); // backscattering amplitude S1 at 180 degrees.
fortran_complex SFORW;
SFORW.r = 0.42;
SFORW.i = 0.42;
fortran_complex SBACK;
SBACK.r = 0.42;
SBACK.i = 0.42;

fortran_complex* TFORW = new fortran_complex[2];
TFORW[0].r = 0.0;
TFORW[0].i = 0.0;
TFORW[1].r = 0.0;
TFORW[1].i = 0.0;
fortran_complex* TBACK = new fortran_complex[2];
TBACK[0].r = 1.0;
TBACK[0].i = 1.0;
TBACK[1].r = 0.0;
TBACK[1].i = 0.0;

float SPIKE;
float* PMOM = new float[4*4] ; // not used.


// Call Wiscombe's mie function to calculate S1 and S2.
MIEV0(&XX, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, 
	  &QEXT, &QSCA, &GQSC, PMOM, &SFORW, &SBACK,  S1, S2, TFORW, TBACK, &SPIKE);


for (int i = 0; i < nangles; i ++) {
    s11[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) + complx_abs(S1[i])*complx_abs(S1[i]));
    s12[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) - complx_abs(S1[i])*complx_abs(S1[i]));
    fortran_complex intermediate = complx_mult(complx_conj(S1[i]), S2[i]);
    s33[i] = intermediate.r;
    s43[i] = intermediate.i;
}

    printf("Beginning simulation\n");
    // Begin simulation
    for (int i = 0; i < nphotons; i ++ ) {
        //printf("Launching photon %d\n", i);
        photon a;
        while (a.alive()) {
            a.move();
            a.drop();
            if (!a.alive()) { break; }
            a.rejection();
            a.scatter();
        }
    }


    printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_R/(nphotons),Q_R/(nphotons),U_R/(nphotons),V_R/(nphotons));	
    printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_T/(nphotons),Q_T/(nphotons),U_T/(nphotons),V_T/(nphotons));	

    delete [] TFORW;
    delete [] TBACK;
	delete [] XMU;
	delete [] PRNT;	
	delete [] S1;
	delete [] S2;
}
