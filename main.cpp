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
double I_R = 0;
double Q_R = 0;	
double U_R = 0;	
double V_R = 0;	

// Stokes vector collector for transmitted photons
double I_T = 0;	
double Q_T = 0;	
double U_T = 0;	
double V_T = 0;	

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

double radius = (2.03/2); // Radius of spherical scatterer 
double wavelength = 0.6328; // wavelength of incident beam 
double rho = 1.152e-4; // Density of spherical scatterers in medium
int nangles = 1000;
double mu_a = 0.0;
double mu_s, slabdepth, albedo; // values to be calculated after Mie

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
float XX= (float) (2 * M_PI * radius / wavelength); // Size parameter
//float XX= 10.0; // Test value
//Complex CREFIN((1.59/1.33),0.0); // fortran_complex refractive index
fortran_complex CREFIN;
CREFIN.r = (1.59/1.33);
CREFIN.i = 0.0;
int PERFCT = 0; // false, refractive index is not infinite.
float MIMCUT = 1e-6; // TODO: is this correct for cutoff of imaginary index, since it's already zero?
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
// OUTPUT
float QEXT, QSCA, GQSC; // Extinction efficiency factor, scattering efficiency factor, and asymmetry factor
// Mie scattering amplitudes (what we need!)
fortran_complex* S1 = new fortran_complex[nangles];
fortran_complex* S2 = new fortran_complex[nangles];
//fortran_complex SFORW(0,0); // Forward scattering amplitude S1 at 0 degrees.
//fortran_complex SBACK(0,0); // backscattering amplitude S1 at 180 degrees.
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

printf("Mie properties: \ndiameter=%5.5f\nmu_s=%5.5f\nrho=%5.5f\nslabdepth=%5.5f\nQSCA=%5.5f\n", radius*2, mu_s, rho, slabdepth, QSCA);

for (int i = 0; i < nangles; i ++) {
    //printf("S1[%d]: %5.5f   S2[%d]: %5.5f\n",i,S1[i].r,i,S2[i].r);
    s11[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) + complx_abs(S1[i])*complx_abs(S1[i]));
    s12[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) - complx_abs(S1[i])*complx_abs(S1[i]));
    fortran_complex intermediate = complx_mult(complx_conj(S1[i]), S2[i]);
    s33[i] = intermediate.r;
    s43[i] = intermediate.i;
}
    printf("Beginning simulation... ");
    fflush(stdout);
    // Begin simulation
    for (int i = 0; i < nphotons; i ++ ) {
        //printf("Launching photon %d\n", i);
        photon a;
        int j = 0;
        while (a.alive()) {
            a.move();
            a.drop();
            if (!a.alive()) { 
                //printf("%5.5f %5.5f %5.5f %5.5f \n",a.S.I,a.S.Q, a.S.U, a.S.V);
               //printf("Photon %d survived %d rounds.\n", i, j);
                break;
            }
            a.rejection();
            a.scatter();
            j ++;
        }
    }

    printf("Simulation done!\n");
    printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_R,Q_R,U_R,V_R);	
    printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_R/(nphotons),Q_R/(nphotons),U_R/(nphotons),V_R/(nphotons));	
    printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_T/(nphotons),Q_T/(nphotons),U_T/(nphotons),V_T/(nphotons));	

    delete [] TFORW;
    delete [] TBACK;
	delete [] XMU;
	delete [] PRNT;	
	delete [] S1;
	delete [] S2;
}
