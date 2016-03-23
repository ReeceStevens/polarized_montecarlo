/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double nphotons = 1e6;
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

double radius = (2.0/2); // Radius of spherical scatterer 
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

double start = clock();

srand(time(NULL));

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
double g = GQSC / QSCA;
printf("\n\n\n***Starting Simulation***\n");
printf("Mie properties: \ndiameter=%5.5f\nmu_s=%5.5f\nrho=%5.5f\nslabdepth=%5.5f\nasymmetry factor=%5.5f\n", radius*2, mu_s, rho, slabdepth, g);

mu_s *= (1-g); // Compensating for backscattering

int grid_res = 100;
const double hw = 7/mu_s;

// Allocate map space
// Profile 1
double** I_Ref = new double*[grid_res];
for (int i = 0; i < grid_res; i ++ ) {
	I_Ref[i] = new double[grid_res];
	for (int j = 0; j < grid_res; j++){
		I_Ref[i][j] = 0;
	}
}

double** Q_Ref = new double*[grid_res];
for (int i = 0; i < grid_res; i ++ ) {
	Q_Ref[i] = new double[grid_res];
	for (int j = 0; j < grid_res; j++){
		Q_Ref[i][j] = 0;
	}
}

double** U_Ref = new double*[grid_res];
for (int i = 0; i < grid_res; i ++ ) {
	U_Ref[i] = new double[grid_res];
	for (int j = 0; j < grid_res; j++){
		U_Ref[i][j] = 0;
	}
}

double** V_Ref = new double*[grid_res];
for (int i = 0; i < grid_res; i ++ ) {
	V_Ref[i] = new double[grid_res];
	for (int j = 0; j < grid_res; j++){
		V_Ref[i][j] = 0;
	}
}

for (int i = 0; i < nangles; i ++) {
    //printf("S1[%d]: %5.5f   S2[%d]: %5.5f\n",i,S1[i].r,i,S2[i].r);
    s11[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) + complx_abs(S1[i])*complx_abs(S1[i]));
    s12[i] = 0.5*(complx_abs(S2[i])*complx_abs(S2[i]) - complx_abs(S1[i])*complx_abs(S1[i]));
    fortran_complex intermediate = complx_mult(complx_conj(S1[i]), S2[i]);
    s33[i] = intermediate.r;
    s43[i] = intermediate.i;
}
    printf("Beginning simulation... \n");
    fflush(stdout);
    // Begin simulation (only doing one orientation for testing purposes)
	for (int j = 0; j < 4; j ++) {
	int num_ref = 0;
    for (int i = 0; i < nphotons; i ++ ) {
        //printf("Launching photon %d\n", i);
        photon a;
		switch(j) {
			case 0:
				a.setStokes(1,1,0,0);
				break;
			case 1:
				a.setStokes(1,-1,0,0);
				break;
			case 2:
				a.setStokes(1,0,1,0);
				break;
			case 3: 
				a.setStokes(1,0,0,1);
				break;
		}
		int k = 0;
        while (a.alive()) {
            a.move();
            a.drop();
            if (!a.alive()) { 
                break;
            }
            a.rejection();
            a.scatter();
            k ++;
			//printf("Step %d\n",k);
        }
		// Move all positions to positive values
		if (a.z < 0) {		
		num_ref ++;
		a.x += hw;
		a.y += hw;
		int idx_x = 0;
		int idx_y = 0;
		// Outside x boundary to the left
		if (a.x < 0) {
			idx_x = 0;
		}
		// Outside x boundary to the right
		else if (a.x > 2*hw) {
			idx_x = grid_res - 1;
		}
		else {
			double dist = a.x / (2*hw);
			dist *= grid_res;
			idx_x = (int) dist;
		}
		// Outside x boundary to the left
		if (a.y < 0) {
			idx_y = 0;
		}
		// Outside x boundary to the right
		else if (a.y > 2*hw) {
			idx_y = grid_res - 1;
		}
		else {
			double dist = a.y / (2*hw);
			dist *= grid_res;
			idx_y = (int) dist;
		}
		I_Ref[idx_x][idx_y] = a.S.I;
		Q_Ref[idx_x][idx_y] = a.S.Q;
		U_Ref[idx_x][idx_y] = a.S.U;
		V_Ref[idx_x][idx_y] = a.S.V;
		}
    }
	printf("Round %d:\n", j);
    printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_R/(nphotons),Q_R/(nphotons),U_R/(nphotons),V_R/(nphotons));	
    printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_T/(nphotons),Q_T/(nphotons),U_T/(nphotons),V_T/(nphotons));	
	printf("Number of reflected photons: %d\n", num_ref);
	I_R = 0;
	Q_R = 0;
	U_R = 0;
	V_R = 0;
	I_T = 0;
	Q_T = 0;
	U_T = 0;
	V_T = 0;
	}

    printf("Simulation done!\n");

	/* Printing Output Matrices */
	FILE* op_matrix;
	op_matrix = fopen("output.m", "w");
	// Output Profile 1
	fprintf(op_matrix, "I_R_1 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", I_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "Q_R_1 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", Q_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "U_R_1 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", U_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "V_R_1 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", V_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");
	// Add MATLAB processing script for easy visualization
	fprintf(op_matrix, "figure();\nsubplot(2,2,1);\nimshow(I_R_1);\ntitle('I_1'); \nsubplot(2,2,2);\nimshow(Q_R_1);\ntitle('Q_1'); \nsubplot(2,2,3);\nimshow(U_R_1);\ntitle('U_1'); \nsubplot(2,2,4);\nimshow(V_R_1);\ntitle('V_1'); \n");

	// Output Profile 2

	fprintf(op_matrix, "I_R_2 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", I_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "Q_R_2 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", Q_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "U_R_2 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", U_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "V_R_2 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", V_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");
	// Add MATLAB processing script for easy visualization
	fprintf(op_matrix, "figure();\nsubplot(2,2,1);\nimshow(I_R_2);\ntitle('I_2'); \nsubplot(2,2,2);\nimshow(Q_R_2);\ntitle('Q_2'); \nsubplot(2,2,3);\nimshow(U_R_2);\ntitle('U_2'); \nsubplot(2,2,4);\nimshow(V_R_2);\ntitle('V_2'); \n");

	// Output Profile 3
	
	fprintf(op_matrix, "I_R_3 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", I_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "Q_R_3 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", Q_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "U_R_3 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", U_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "V_R_3 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", V_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");
	// Add MATLAB processing script for easy visualization
	fprintf(op_matrix, "figure();\nsubplot(2,2,1);\nimshow(I_R_3);\ntitle('I_3'); \nsubplot(2,2,2);\nimshow(Q_R_3);\ntitle('Q_3'); \nsubplot(2,2,3);\nimshow(U_R_3);\ntitle('U_3'); \nsubplot(2,2,4);\nimshow(V_R_3);\ntitle('V_3'); \n");

	// Output Profile 4

	fprintf(op_matrix, "I_R_4 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", I_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "Q_R_4 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", Q_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "U_R_4 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", U_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");

	fprintf(op_matrix, "V_R_4 = [");
	for (int i = 0; i < grid_res; i ++) {
		for (int j = 0; j < grid_res; j ++) {
			fprintf(op_matrix, "%f ", V_Ref[i][j]);	
		}
		fprintf(op_matrix,";\n");
	}
	fprintf(op_matrix,"];\n");
	// Add MATLAB processing script for easy visualization
	fprintf(op_matrix, "figure();\nsubplot(2,2,1);\nimshow(I_R_4);\ntitle('I_4'); \nsubplot(2,2,2);\nimshow(Q_R_4);\ntitle('Q_4'); \nsubplot(2,2,3);\nimshow(U_R_4);\ntitle('U_4'); \nsubplot(2,2,4);\nimshow(V_R_4);\ntitle('V_4'); \n");

	fclose(op_matrix);

	double finish = clock();

	printf("Total time: %5f seconds\n", (finish-start)/CLOCKS_PER_SEC);
	printf("\n\n");
    delete [] TFORW;
    delete [] TBACK;
	delete [] XMU;
	delete [] PRNT;	
	delete [] S1;
	delete [] S2;
}
