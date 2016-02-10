/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
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

//struct { float r,i; } fortran_complex;
//typedef struct fortran_complex fortran_complex;

float pi = 3.1415926536;
double radius = 1.00; // TODO: check units of radius. I think microns?
double wavelength = 0.600; // TODO: microns?
int nangles = 1000;

extern "C" {
void MIEV0_(float*, Complex*, int*, float*, int*, float*, float*, int*, int*, int*, int*, float*, 
			float*, float*, Complex*, Complex*, Complex*, Complex*, Complex*, Complex*, 
			float*, float**);
}


int main() {


// Mie scattering variables. Passed into Wiscombe's mie scattering subroutine.
// Initial values taken from Ramella et. al. model.
// INPUT
float XX= (float) (2 * pi * radius / wavelength); // Size parameter
Complex CREFIN((1.59/1.33),0.0); // Complex refractive index
int PERFCT = 0; // false, refractive index is not infinite.
float MIMCUT = 0; // TODO: is this correct for cutoff of imaginary index, since it's already zero?
int ANYANG = 0; // Angles are monotone increasing and mirror symmetric around 90 degrees
float NUMANG = (float) nangles; // number of angles to be evaluated for S1 and S2.
float* XMU = new float[nangles];
// fill out angle array
for(int i = 0; i <nangles; i ++) {
	XMU[i] = pi * i / nangles;
}
int NMOM = 0; // Avoid calculation of PMOM
int IPOLZN = 0; // Not used.
int MOMDIM = 0; // Not used.
int* PRNT = new int[2];
PRNT[0] = 0;
PRNT[1] = 0;
// OUTPUT
float QEXT, QSCA, GQSC; // Extinction efficiency factor, scattering efficiency factor, and asymmetry factor
// Mie scattering amplitudes (what we need!)
Complex* S1 = new Complex[nangles];
Complex* S2 = new Complex[nangles];
Complex SFORW(0,0); // Forward scattering amplitude S1 at 0 degrees.
Complex SBACK(0,0); // backscattering amplitude S1 at 180 degrees.
Complex* TFORW = new Complex[2];
Complex* TBACK = new Complex[2];
float SPIKE;
float** PMOM; // not used.


// Call Wiscombe's mie function to calculate S1 and S2.
MIEV0_(&XX, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, 
	  &QEXT, &QSCA, &GQSC, S1, S2, &SFORW, &SBACK, TFORW, TBACK, &SPIKE, PMOM);



	// Check Mie by printing
	for (int i = 0; i < NUMANG; i ++) {
		printf("S1[%d]: %f, %f. S2[%d]: %f, %f \n", i, S1[i].Re(), S1[i].Im(), i, S2[i].Re(), S2[i].Im());
	}


	delete [] XMU;
	delete [] PRNT;	
	delete [] S1;
	delete [] S2;
}
