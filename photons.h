#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "Vector.h"
#include "Complex.h"
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
//
// Validation using newest Ramella paper. compare to ramella output for validation.


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
	stokes_v(const stokes_v& that) {
		this->I = that.I;
		this->Q = that.Q;
		this->U = that.U;
		this->V = that.V;
	}
    void scalar_mult(double x) { I = x*I; Q = x*Q; U = x*U; V = x*V;}
	void rotate_stokes(double phi) {
		double cos_2phi = cos(2*phi);
		double sin_2phi = sin(2*phi);
		double old_q = this->Q;
		double old_u = this->U;
		this->Q = old_q*cos_2phi + old_u*sin_2phi;
		this->U = -old_q*sin_2phi + old_u*cos_2phi;
	}
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
	double alpha = 0; // scattering angles
	double beta = 0; // scattering angles
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
     * or out of bounds. Adjust weights accordingly. If photon
	 * is dead, re-orient it to detector and record stokes vector.
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
		if (z <= 0) { // Photon is reflected
			// Rotate axes back to reference frame
			double x_on_z = U.i*V.j - U.j*V.i;
			// Empty case: the vector is already correctly oriented.
			if !((abs(x_on_z) < 1e-8) && (abs(V.k) < 1e-8)) {
				double phi = atan2(V.k,-x_on_z);
				S.rotate_stokes(phi);
				phi = atan2(U.j, -U.i);
				S.rotate_stokes(phi);
			}
			I_R += S.I;
			Q_R += S.Q;
			U_R += S.U;
			V_R += S.V;
			state = DEAD;
			
			// TODO: Quantize x and y location, place stokes values in matrix
					
		}
		else if (z >= SLABSIZE_Z) { // Photon is transmitted
			// Rotate axes back to reference frame
			double x_on_z = U.i*V.j - U.j*V.i;
			// Empty case: the vector is already correctly oriented.
			if !((abs(x_on_z) < 1e-8) && (abs(V.k) < 1e-8)) {
				double phi = atan2(V.k,-x_on_z);
				S.rotate_stokes(phi);
				phi = atan2(U.j, U.i);
				S.rotate_stokes(phi);
			}
			I_T += S.I;
			Q_T += S.Q;
			U_T += S.U;
			V_T += S.V;
			state = DEAD;

			// TODO: Quantize x and y location, place stokes values in matrix
		}
    }

	void rejection(void) {
		alpha = 1.00;
		beta = 3.14;
		return;
		/*double* angles = new double[2];
		angles[0] = 1.00;
		angles[1] = 3.14;
		return angles;*/
	}

	void rotate_about_vector(Vector& trajectory, Vector& axis, double angle) {
		int columns = 3;
		int rows = 3;	
		double** r_euler = new double* [columns];
		for (int i = 0; i < columns; i += 1) {
			r_euler[i] = new double [rows];
		}

		// Define a rotational matrix R_{euler}
		r_euler[0][0] = axis.i*axis.i*(1- cos(angle));
		r_euler[1][0] = axis.j*axis.i*(1- cos(angle)) - axis.k*sin(angle);	
		r_euler[2][0] = axis.k*axis.i*(1- cos(angle)) + axis.j*sin(angle);
		r_euler[0][1] = axis.i*axis.j*(1- cos(angle)) + axis.k*sin(angle);	
		r_euler[1][1] = axis.j*axis.j*(1- cos(angle)) + cos(angle);	
		r_euler[2][1] = axis.j*axis.k*(1- cos(angle)) - axis.i*sin(angle);	
		r_euler[0][2] = axis.i*axis.k*(1- cos(angle)) - axis.j*sin(angle);	
		r_euler[1][2] = axis.j*axis.k*(1- cos(angle)) + axis.i*sin(angle);	
		r_euler[2][2] = axis.k*axis.k*(1- cos(angle)) + cos(angle);	

		// Make vectors into arrays for easier multiplication
		double* v = new double [3];	
		v[0] = trajectory.i;	
		v[1] = trajectory.j;	
		v[2] = trajectory.k;	

		double* new_v = new double [3];	
		new_v[0] = 0;	
		new_v[1] = 0;	
		new_v[2] = 0;

		// Multiply trajectory by rotational matrix.
		for (int i = 0; i < 3; i += 1) {
			for (int j = 0; j < 3; j += 1) {
				new_v[i] += v[j] * r_euler[j][i];
			}			
		}

		// Update trajectory vector
		trajectory.i = new_v[0];
		trajectory.j = new_v[1];
		trajectory.k = new_v[2];	

		delete [] r_euler;
		delete [] v;
	}

    /*
     * TODO:
     * scatter(void) - 
     * Performs scattering based on the polarizing angle and scattering
     * particles in the media. For spheres, use Mie scattering theory.
     */
    void scatter(void) {
		// Determine scattering angles by rejection method.
		/*double* scat_angles = rejection();
		double alpha = scat_angles[0];
		double beta = scat_angles[1];*/

		// Rotate V about U by alpha
		rotate_about_vector(this.V, this.U, alpha);

		// Rotate U about new V by beta
		rotate_about_vector(this.U, this.V, beta);

		// TODO: Adjust the Stokes vector
					


    }
};

#endif
