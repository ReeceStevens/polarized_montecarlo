#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "Vector.h"
#include "mie.h"

#ifndef __photons__h
#define __photons__h 1

#define ALIVE 1
#define DEAD 0

// TODO: 
// make the effects of diattenuation optional (turn on or off)
// decide whether or not to do the sphere or cylinder at every scattering encounter
//
// Get spherical model working using ramella code or other lit
// test decision making by choosing two different diameters of spherical scatterers
//
// Make modular? more general structure, additional components for later montecarlo
//
// swappable sizes, polarization or not, birefringence or not, etc.
//
// new phase functions or not


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

/* Properties of the light */


/* Properties of our model (taken from Ramella et. al., 2005) */
#define WEIGHT_CUTOFF 0.01 // From Ramella et. al.
const double survival_multiplier = 10;
const double survival_chance = 1 / survival_multiplier;

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
    void scalar_mult(double x) { I = x*I; Q = x*Q; U = x*U; V = x*V;}
};

struct reference_frame {
private:
	Vector V;
	Vector U;
public:	
	Vector returnThird() {
		return Vector.cross_prod(V,U);
	}
};

struct photon {
private:
    stokes_v S;
    Vector V; // Orientation vector
    Vector U; // Direction cosine vector
    int32_t x;
    int32_t y;
    int32_t z;
    uint32_t weight;
    bool state; // Alive or dead
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
        state = ALIVE;
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
        // Check weight
        if (weight < WEIGHT_CUTOFF) {  
            // Russian Roulette for photon absorption
            double rand = rand_num();
            if (rand > survival_chance) {
                weight = 0;  
                state = DEAD;
            } else {    
                weight *= survival_multiplier; 
            } 
        }
        // Check boundaries
        if ((x > SLABSIZE_X) || (x < 0) || (y > SLABSIZE_Y) || (y < 0) || (z > SLABSIZE_Z) || (z < 0)) {
            S.scalar_mult(weight);          
            state = DEAD;
        }
    }

	double* rejection(void) {
		double* angles = new double[2];
		angles[0] = 1.00;
		angles[1] = 3.14;
		return angles;
	}

    /*
     * TODO:
     * scatter(void) - 
     * Performs scattering based on the polarizing angle and scattering
     * particles in the media. For spheres, use Mie scattering theory.
     */
    void scatter(void) {
		// Determine scattering angles by rejection method.
		double* scat_angles = rejection();
		double alpha = scat_angles[0];
		double beta = scat_angles[1];
		// Rotate V about U
		int columns = 3;
		int rows = 3;	
		double** data = new double* [columns];
		for (int i = 0; i < columns; i += 1) {
			data[i] = new double [rows];
		}
		// Define a rotational matrix R_{euler}
		data[0][0] = U.i*U.i*(1- cos(alpha));
		data[1][0] = U.j*U.i*(1- cos(alpha)) - U.k*sin(alpha);	
		data[2][0] = U.k*U.i*(1- cos(alpha)) + U.j*sin(alpha);
		data[0][1] = U.i*U.j*(1- cos(alpha)) + U.k*sin(alpha);	
		data[1][1] = U.j*U.j*(1- cos(alpha)) + cos(alpha);	
		data[2][1] = U.j*U.k*(1- cos(alpha)) - U.i*sin(alpha);	
		data[0][2] = U.i*U.k*(1- cos(alpha)) - U.j*sin(alpha);	
		data[1][2] = U.j*U.k*(1- cos(alpha)) + U.i*sin(alpha);	
		data[2][2] = U.k*U.k*(1- cos(alpha)) + cos(alpha);	

		Vector newV = new Vector	

		// TODO: does this delete the individual mallocs inside the matrix?
		delete [] scat_angles;
		delete [] data;
    }
};

#endif
