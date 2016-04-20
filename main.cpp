/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Declare variables we're grabbing from config.h
double I_R, Q_R, U_R, V_R, I_T, Q_T, U_T, V_T;
double mu_a, mu_s, slabdepth, albedo, g;
int nangles, nphotons;
double *s11, *s12, *s33, *s43;
double **I_Ref_H, **I_Ref_V, **I_Ref_P, **I_Ref_M, **I_Ref_R, **I_Ref_L;
double **Q_Ref_H, **Q_Ref_V, **Q_Ref_P, **Q_Ref_M, **Q_Ref_R, **Q_Ref_L;
double **U_Ref_H, **U_Ref_V, **U_Ref_P, **U_Ref_M, **U_Ref_R, **U_Ref_L;
double **V_Ref_H, **V_Ref_V, **V_Ref_P, **V_Ref_M, **V_Ref_R, **V_Ref_L;

#include "config.h"
#include "photons.h"
#include "printout.h"


int main() {

	config();
	double start = clock();

	srand(time(NULL));

	printf("\n\n\n***Starting Simulation***\n");
	//printf("Mie properties: \ndiameter=%5.5f\nmu_s=%5.5f\nrho=%5.5f\nslabdepth=%5.5f\nasymmetry factor=%5.5f\n", radius*2, mu_s, rho, slabdepth, g);
	printf("Mie properties: \nmu_s=%5.5f\nslabdepth=%5.5f\ng=%5.5f\n", mu_s, slabdepth,g);

	int grid_res = 100;
	const double hw = 7/(mu_s);

	printf("Beginning simulation... \n");
	fflush(stdout);
	// Begin simulation (only doing one orientation for testing purposes)
	for (int j = 0; j < 6; j ++) {
		int num_ref = 0;
		for (int i = 0; i < nphotons; i ++ ) {
			if (i%100000 == 0) { printf("\r%d percent complete!",(i/10000)); }
			fflush(stdout);
			photon a = photon();
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
					a.setStokes(1,0,-1,0);
					break;
				case 4:
					a.setStokes(1,0,0,1);
					break;
				case 5:
					a.setStokes(1,0,0,-1);
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
				a.quaternion_scatter();
				k ++;
			}
			// Move all positions to positive values
			if (a.z < 0) {		
			num_ref ++;
			int idx_x = 0;
			int idx_y = 0;
			// Outside x boundary to the left
			if (a.x < -hw) {
				idx_x = 0;
			}
			// Outside x boundary to the right
			else if (a.x > hw) {
				idx_x = grid_res - 1;
			}
			// TODO: What's the justification for this out-of-bounds behavior?
			else {
				double dist = fabs(a.x + hw) * grid_res / (2*hw);
				idx_x = (int) dist;
			}
			// Outside x boundary to the left
			if (a.y < -hw) {
				idx_y = 0;
			}
			// Outside x boundary to the right
			else if (a.y > hw) {
				idx_y = grid_res - 1;
			}
			else {
				double dist = fabs(a.y + hw) * grid_res / (2*hw);
				idx_y = (int) dist;
			}
			//a.S.normalize();
			switch(j) {
				case 0:
					I_Ref_H[idx_y][idx_x] += a.S.I;
					Q_Ref_H[idx_y][idx_x] += a.S.Q;
					U_Ref_H[idx_y][idx_x] += a.S.U;
					V_Ref_H[idx_y][idx_x] += a.S.V;
					break;
				case 1:
					I_Ref_V[idx_y][idx_x] += a.S.I;
					Q_Ref_V[idx_y][idx_x] += a.S.Q;
					U_Ref_V[idx_y][idx_x] += a.S.U;
					V_Ref_V[idx_y][idx_x] += a.S.V;
					break;
				case 2:
					I_Ref_P[idx_y][idx_x] += a.S.I;
					Q_Ref_P[idx_y][idx_x] += a.S.Q;
					U_Ref_P[idx_y][idx_x] += a.S.U;
					V_Ref_P[idx_y][idx_x] += a.S.V;
					break;
				case 3:
					I_Ref_M[idx_y][idx_x] += a.S.I;
					Q_Ref_M[idx_y][idx_x] += a.S.Q;
					U_Ref_M[idx_y][idx_x] += a.S.U;
					V_Ref_M[idx_y][idx_x] += a.S.V;
					break;
				case 4:
					I_Ref_R[idx_y][idx_x] += a.S.I;
					Q_Ref_R[idx_y][idx_x] += a.S.Q;
					U_Ref_R[idx_y][idx_x] += a.S.U;
					V_Ref_R[idx_y][idx_x] += a.S.V;
					break;
				case 5:
					I_Ref_L[idx_y][idx_x] += a.S.I;
					Q_Ref_L[idx_y][idx_x] += a.S.Q;
					U_Ref_L[idx_y][idx_x] += a.S.U;
					V_Ref_L[idx_y][idx_x] += a.S.V;
					break;
			}
			}
    }
	printf("\nRound %d:\n", j);
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
	printout();
	double finish = clock();
	printf("Total time: %5f seconds\n", (finish-start)/CLOCKS_PER_SEC);
	printf("\n\n");

}
