#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "Vector.h"

#ifndef __photons__h
#define __photons__h 1

/* Properties of the medium */
#define SLABSIZE_X 10
#define SLABSIZE_Y 10
#define SLABSIZE_Z 30
#define MU_S 0.4
#define MU_A 0.5
#define ALBEDO ((MU_S) / (MU_S + MU_A))
#define MU_TOTAL (MU_A + MU_S)
const double mu_total = MU_A + MU_S;
const double albedo = MU_S / (mu_total);

/* Properties of our model */
#define WEIGHT_CUTOFF 0.01 // From Ramella et. al


/* Some random number generator that produces a uniform 
 * distribution of numbers from (0,1].
 */
double rand_num();

/* Track polarization state of photon via Stokes vector */
struct stokes_v {
private:
    double I;
    double Q;
    double U;
    double V;
public:
    stokes_v(void) { I = 0; Q = 0; U = 0; V = 0; }
    stokes_v(double I, double Q, double U, double V) : I(I), Q(Q), U(U), V(V) { }; 
};

struct photon {
private:
    stokes_v S;
    Vector V; // Orientation vector
    Vector U; // Direction cosine vector
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t weight;
public:
    // Default constructor.
    photon(void) {
        weight = 1;
        // Start at center and top of slab
        x = SLABSIZE_X / 2;
        y = SLABSIZE_Y / 2;
        z = 0;
        Vector V(0,0,1);
        Vector U(0,0,1);
    } 

    /* 
     * move(void) -
     * Adjust position of photon in space based upon
     * orientation and properties of the medium.
     */
    void move(void) {
        double deltaS = -log(rand_num()) / mu_total;
        x = x + U.i * deltaS;
        y = y + U.j * deltaS;
        z = z + U.k * deltaS;
    }

    /* 
     * drop(void) -
     * Drop weight from the photon and determine if it is dead
     * or out of bounds. Adjust weights accordingly.
     */
    void drop(void) {
        weight *= albedo; 
        if (weight < WEIGHT_CUTOFF) { roulette();}       
    }
};

#endif
