#ifndef __mie_h__
#define __mie_h__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* globals */
double I_R, Q_R, U_R, V_R, I_T, Q_T, U_T, V_T;
double mu_a, mu_s, slabdepth, albedo, g;
double *s11, *s12, *s33, *s43;
int nangles, nphotons, grid_res;


struct fortran_complex{ float r,i; };

/* Mie Setup */
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

extern "C" void MIEV0(float* XX, fortran_complex* CREFIN, int* PERFCT, float* MIMCUT, int* ANYANG, 
                    int* NUMANG, 
                   float* XMU, int* NMOM, int* IPOLZN, int* MOMDIM, int* PRNT, float* QEXT, 
			       float* QSCA, float* GQSC, float* PMOM, fortran_complex* SFORW, fortran_complex* SBACK,
                   fortran_complex* S1, fortran_complex* S2, 
                   fortran_complex* TFORW, fortran_complex* TBACK, float* SPIKE );

void callWiscombeMie() {
	nphotons = 1e6;
	/* Create Global Variables */
	// Stokes vector collector for reflected photons
	I_R = 0;
	Q_R = 0;
	U_R = 0;
	V_R = 0;

	// Stokes vector collector for transmitted photons
	I_T = 0;
	Q_T = 0;
	U_T = 0;
	V_T = 0;

	//double radius = (2.0/2); // Radius of spherical scatterer 
	double radius = (2.02/2); // Radius of spherical scatterer 
	//double wavelength = 0.6328; // wavelength of incident beam 
	double wavelength = 0.543; // wavelength of incident beam 
	double rho = 1.152e-4; // Density of spherical scatterers in medium
	nangles = 1000;
	mu_a = 0.0;

	s11 = new double[nangles];
	s12 = new double[nangles];
	s33 = new double[nangles];
	s43 = new double[nangles];

	// Mie scattering variables. Passed into Wiscombe's mie scattering subroutine.
	// Initial values taken from Ramella et. al. model.
	// INPUT
	fortran_complex CREFIN;
	CREFIN.r = (1.59/1.33);
	CREFIN.i = 0.0;
	int PERFCT = 0; // false, refractive index is not infinite.
	float MIMCUT = 1e-6; 
	float XX= (float) (2 * M_PI * radius / (wavelength/1.33)); // Size parameter
	//float XX= 10.0; // Test value
	int ANYANG = 0; // Angles are monotone increasing and mirror symmetric around 90 degrees
	int NUMANG = nangles; // number of angles to be evaluated for S1 and S2.
	float* XMU = new float[nangles];
	// fill out angle array
	for(int i = 0; i <nangles; i ++) {
		XMU[i] = (float) cos(M_PI * i / nangles);
	}
	int NMOM = 0; // Avoid calculation of PMOM
	int IPOLZN = -1; // Not used.
	int MOMDIM = 1; // Not used.
	int* PRNT = new int[2];
	PRNT[0] = 0;
	PRNT[1] = 0;
	/* OUTPUT */
	// Extinction efficiency factor, scattering efficiency factor, and asymmetry factor
	float QEXT, QSCA, GQSC; 
	// Mie scattering amplitudes (what we need!)
	fortran_complex* S1 = new fortran_complex[nangles];
	fortran_complex* S2 = new fortran_complex[nangles];
	fortran_complex SFORW;
	SFORW.r = 0.0;
	SFORW.i = 0.0;
	fortran_complex SBACK;
	SBACK.r = 0.0;
	SBACK.i = 0.0;

	fortran_complex* TFORW = new fortran_complex[2];
	TFORW[0].r = 0.0;
	TFORW[0].i = 0.0;
	TFORW[1].r = 0.0;
	TFORW[1].i = 0.0;
	fortran_complex* TBACK = new fortran_complex[2];
	TBACK[0].r = 0.0;
	TBACK[0].i = 0.0;
	TBACK[1].r = 0.0;
	TBACK[1].i = 0.0;

	float SPIKE;
	float* PMOM = new float[4*4] ; // not used.

	// Call Wiscombe's mie function to calculate S1 and S2.
	MIEV0(&XX, &CREFIN, &PERFCT, &MIMCUT, &ANYANG, &NUMANG, XMU, &NMOM, &IPOLZN, &MOMDIM, PRNT, 
	  &QEXT, &QSCA, &GQSC, PMOM, &SFORW, &SBACK,  S1, S2, TFORW, TBACK, &SPIKE);

	/* Properties of the medium */
	mu_s = QSCA*M_PI*radius*radius*rho*1e4; // inverse cm
	slabdepth = 4/mu_s;
	albedo = mu_s / (mu_s + mu_a);
	g = GQSC / QSCA;
	grid_res = 100;

	for (int i = 0; i < nangles; i ++) {
		s11[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) + complx_abs(S1[i])*complx_abs(S1[i]));
		s12[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) - complx_abs(S1[i])*complx_abs(S1[i]));
		fortran_complex intermediate = complx_mult(complx_conj(S1[i]), S2[i]);
		s33[i] = intermediate.r;
		s43[i] = intermediate.i;
	}
    delete [] TFORW;
    delete [] TBACK;
	delete [] XMU;
	delete [] PRNT;
    delete [] PMOM;
    delete [] S1;
    delete [] S2;
}
#endif
