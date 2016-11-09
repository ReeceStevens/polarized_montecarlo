/*
 * config.h
 * Parse user configuration files to produce variables for Monte Carlo.
 * Currently, just holds the setup for testing. Future work: expand this
 * file to parse JSON for config parameters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mie.h"

// These variables are coming from main 
extern scat_props_t* mie_props;
extern double I_R, Q_R, U_R, V_R, I_T, Q_T, U_T, V_T;
extern double mu_a, mu_s, slabdepth, albedo, g;
extern int nphotons, nangles, grid_res;
extern double *s11, *s12, *s33, *s43;
extern double **I_Ref_H, **I_Ref_V, **I_Ref_P, **I_Ref_M, **I_Ref_R, **I_Ref_L;
extern double **Q_Ref_H, **Q_Ref_V, **Q_Ref_P, **Q_Ref_M, **Q_Ref_R, **Q_Ref_L;
extern double **U_Ref_H, **U_Ref_V, **U_Ref_P, **U_Ref_M, **U_Ref_R, **U_Ref_L;
extern double **V_Ref_H, **V_Ref_V, **V_Ref_P, **V_Ref_M, **V_Ref_R, **V_Ref_L;

static void getMapSpace(void);

int config(void) {
    nphotons = 1e6;

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

    grid_res = 100;

    mie_props = mie_setup();
    getMapSpace();

    return 0;
}

static void getMapSpace() {
    // Profile 1
    I_Ref_H = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_H[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_H[i][j] = 0;
        }
    }

    Q_Ref_H = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_H[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_H[i][j] = 0;
        }
    }

    U_Ref_H = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_H[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_H[i][j] = 0;
        }
    }

    V_Ref_H = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_H[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_H[i][j] = 0;
        }
    }

    // Profile 2
    I_Ref_V = new double* [grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_V[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_V[i][j] = 0;
        }
    }

    Q_Ref_V = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_V[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_V[i][j] = 0;
        }
    }

    U_Ref_V = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_V[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_V[i][j] = 0;
        }
    }

    V_Ref_V = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_V[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_V[i][j] = 0;
        }
    }

    // Profile 3
    I_Ref_M = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_M[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_M[i][j] = 0;
        }
    }

    Q_Ref_M = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_M[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_M[i][j] = 0;
        }
    }

    U_Ref_M = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_M[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_M[i][j] = 0;
        }
    }

    V_Ref_M = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_M[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_M[i][j] = 0;
        }
    }

    // Profile 4
    I_Ref_P = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_P[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_P[i][j] = 0;
        }
    }

    Q_Ref_P = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_P[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_P[i][j] = 0;
        }
    }

    U_Ref_P = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_P[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_P[i][j] = 0;
        }
    }

    V_Ref_P = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_P[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_P[i][j] = 0;
        }
    }

    // Profile 5
    I_Ref_R = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_R[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_R[i][j] = 0;
        }
    }

    Q_Ref_R = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_R[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_R[i][j] = 0;
        }
    }

    U_Ref_R = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_R[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_R[i][j] = 0;
        }
    }

    V_Ref_R = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_R[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_R[i][j] = 0;
        }
    }

    // Profile 6
    I_Ref_L = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        I_Ref_L[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            I_Ref_L[i][j] = 0;
        }
    }

    Q_Ref_L = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        Q_Ref_L[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            Q_Ref_L[i][j] = 0;
        }
    }

    U_Ref_L = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        U_Ref_L[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            U_Ref_L[i][j] = 0;
        }
    }

    V_Ref_L = new double*[grid_res];
    for (int i = 0; i < grid_res; i ++ ) {
        V_Ref_L[i] = new double[grid_res];
        for (int j = 0; j < grid_res; j++){
            V_Ref_L[i][j] = 0;
        }
    }
}
