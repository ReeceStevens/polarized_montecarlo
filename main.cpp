/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Complex.h"

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
//typedef struct fortran_complex fortran_complex;

float pi = 3.1415926536;
double radius = 1.00; // TODO: check units of radius. I think microns?
double wavelength = 0.600; // TODO: microns?  
int nangles = 1000;

extern "C" void MIEV0(float* XX, fortran_complex* CREFIN, int* PERFCT, float* MIMCUT, int* ANYANG, 
                    int* NUMANG, 
                   float* XMU, int* NMOM, int* IPOLZN, int* MOMDIM, int* PRNT, float* QEXT, 
			       float* QSCA, float* GQSC, float* PMOM, fortran_complex* SFORW, fortran_complex* SBACK,
                   fortran_complex* S1, fortran_complex* S2, 
                   fortran_complex* TFORW, fortran_complex* TBACK, float* SPIKE );

int main() {


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

    delete [] TFORW;
    delete [] TBACK;
	delete [] XMU;
	delete [] PRNT;	
	delete [] S1;
	delete [] S2;
}
