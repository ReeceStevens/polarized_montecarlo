/* main.cpp
 *
 * Main file for the polarized monte carlo simulation
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "scat.h"

// Declare variables we're grabbing from config.h
scat_props_t* mie_props;
double I_R, Q_R, U_R, V_R, I_T, Q_T, U_T, V_T;
double mu_a, slabdepth, albedo, g;
int nangles, nphotons, grid_res;
double *s11, *s12, *s33, *s43;
double **I_Ref_H, **I_Ref_V, **I_Ref_P, **I_Ref_M, **I_Ref_R, **I_Ref_L;
double **Q_Ref_H, **Q_Ref_V, **Q_Ref_P, **Q_Ref_M, **Q_Ref_R, **Q_Ref_L;
double **U_Ref_H, **U_Ref_V, **U_Ref_P, **U_Ref_M, **U_Ref_R, **U_Ref_L;
double **V_Ref_H, **V_Ref_V, **V_Ref_P, **V_Ref_M, **V_Ref_R, **V_Ref_L;

#include "config.h"
#include "mie.h"
#include "photons.h"
#include "printout.h"

void updateProgressBar(int photonNumber) {
    if ((photonNumber%100000) == 0) { 
        int progress = photonNumber/100000;
        printf("\rProgress: [");
        for (int k = 0; k < progress; k++){
            printf("••••");
        }
        for (int k = 0; k < (10-progress); k++){
            printf("----");
        }
        printf("] ");
        printf("%d%%",progress*10); 
        fflush(stdout);
    }
}

void setStokesVector(photon a, int simulationSet) {
    switch(simulationSet) {
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
}

int launchPhoton(photon a) {
    int k = 0;
    while (a.alive()) {
        a.move(mie_props->mu_s);
        a.drop();
        if (!a.alive()) { 
            break;
        }
        a.rejection();
        a.quaternion_scatter();
        k ++;
    }
    return k;
}

void logPhoton(photon a, const int simulationSet, const double hw) {
    if (a.z < 0) {	// If reflected, mark on surface map
        int idx_x = 0;
        int idx_y = 0;
        // X index
        if (a.x > hw) {
            idx_x = grid_res - 1;
        } else if (a.x > -hw) {
            double dist = fabs(a.x + hw) * (grid_res-1) / (2*hw);
            idx_x = (int) dist;
        }
        // Y index
        if (a.y > hw) {
            idx_y = grid_res - 1;
        } else if (a.y > -hw) {
            double dist = fabs(a.y + hw) * (grid_res-1) / (2*hw);
            idx_y = (int) dist;
        }
        switch(simulationSet) {
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

int main(int argc, char* argv[]) {
    printf("Beginning\n");
    config();
    printf("Configured\n");
    double start = clock();
    printf("Clock started\n");
    srand(time(NULL));
    printf("RNG Primed\n");
    const double avg_ang = 0;
    printf("mu_s function setting\n");
    /* double (*mu_s)(double ang) = mie_props->mu_s; */
    printf("mu_s function set\n");
    const double hw = 7/(mie_props->mu_s(avg_ang));
    printf("hw calculated\n");

    printf("\n\n\n***Starting Simulation***\n");
    printf("Mie properties: \nmu_s=%5.5f\nslabdepth=%5.5f\ng=%5.5f\n", mie_props->mu_s(avg_ang), slabdepth,g);
    printf("Beginning simulation... \n\n");
    for (int j = 0; j < 6; j ++) {
        int num_ref = 0;
        double stotal = 0;
        for (int i = 0; i < nphotons; i ++ ) {
            updateProgressBar(i);
            // Single Photon Simulation
            photon a = photon();
            // Choose incident polarization state
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
            // Life cycle of a photon
            launchPhoton(a);
            // Photon has died or exited tissue. Mark location and tally.
            logPhoton(a, j, hw);
    }
    printf("\rProgress: [••••••••••••••••••••••••••••••••••••••••] 100%%\n");
    printf("\nRound %d:\n", j);
    printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_R/(nphotons),Q_R/(nphotons),U_R/(nphotons),V_R/(nphotons));    
    printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",I_T/(nphotons),Q_T/(nphotons),U_T/(nphotons),V_T/(nphotons));    
    printf("Number of reflected photons: %d\n", num_ref);
    printf("Average total step size: %5.5f\n\n", stotal/nphotons);

    // Reset for next set of photons
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
    // If no output file name given, just use default.
    if (argc < 2) {
        printout("");
    }
    // Otherwise, use given prefix
    else {
        printout(argv[1]);
    }
    double finish = clock();
    printf("Total time: %5f seconds\a\n", (finish-start)/CLOCKS_PER_SEC);
    printf("\n\n");
}
