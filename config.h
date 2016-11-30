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

#define GRID_RES 100

/* globals */
double I_Ref_H[GRID_RES][GRID_RES];
double I_Ref_V[GRID_RES][GRID_RES];
double I_Ref_P[GRID_RES][GRID_RES];
double I_Ref_M[GRID_RES][GRID_RES];
double I_Ref_R[GRID_RES][GRID_RES];
double I_Ref_L[GRID_RES][GRID_RES];

double Q_Ref_H[GRID_RES][GRID_RES];
double Q_Ref_V[GRID_RES][GRID_RES];
double Q_Ref_P[GRID_RES][GRID_RES];
double Q_Ref_M[GRID_RES][GRID_RES];
double Q_Ref_R[GRID_RES][GRID_RES];
double Q_Ref_L[GRID_RES][GRID_RES];

double U_Ref_H[GRID_RES][GRID_RES];
double U_Ref_V[GRID_RES][GRID_RES];
double U_Ref_P[GRID_RES][GRID_RES];
double U_Ref_M[GRID_RES][GRID_RES];
double U_Ref_R[GRID_RES][GRID_RES];
double U_Ref_L[GRID_RES][GRID_RES];

double V_Ref_H[GRID_RES][GRID_RES];
double V_Ref_V[GRID_RES][GRID_RES];
double V_Ref_P[GRID_RES][GRID_RES];
double V_Ref_M[GRID_RES][GRID_RES];
double V_Ref_R[GRID_RES][GRID_RES];
double V_Ref_L[GRID_RES][GRID_RES];

int config(void) {
    callWiscombeMie();
	return 0;
}
