#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "Vector.h"

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
#define SLABSIZE_Z 1
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
double rand_num() {
   double r = ((double) rand() / RAND_MAX); 
   return r;
}

/* Track polarization state of photon via Stokes vector */
struct stokes_v {
public:
    double I;
    double Q;
    double U;
    double V;
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
    void normalize() {
        if (I = 0) { return; }
        Q = Q/I;
        U = U/I;
        V = V/I;
        I = 1.0;
    }
};


struct photon {
private:
    stokes_v S;
    //Vector V; // Orientation vector
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
        //Vector V(0,0,1);
        Vector U(0,0,1);
        state = ALIVE;
    } 

    bool alive() {
        return state; 
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
			//double x_on_z = U.i*V.j - U.j*V.i;
			// Empty case: the vector is already correctly oriented.
			//if !((abs(x_on_z) < 1e-8) && (abs(V.k) < 1e-8)) {
				//double phi = atan2(V.k,-x_on_z);
				//S.rotate_stokes(phi);
				double phi = atan2(U.j, -U.i);
				S.rotate_stokes(phi);
			//}
			I_R += S.I;
			Q_R += S.Q;
			U_R += S.U;
			V_R += S.V;
			state = DEAD;
			
			// TODO: Quantize x and y location, place stokes values in matrix
					
		}
		else if (z >= SLABSIZE_Z) { // Photon is transmitted
			// Rotate axes back to reference frame
			//double x_on_z = U.i*V.j - U.j*V.i;
			// Empty case: the vector is already correctly oriented.
			//if !((abs(x_on_z) < 1e-8) && (abs(V.k) < 1e-8)) {
				//double phi = atan2(V.k,-x_on_z);
				//S.rotate_stokes(phi);
			double phi = -atan2(U.j, U.i);
			S.rotate_stokes(phi);
			//}
			I_T += S.I;
			Q_T += S.Q;
			U_T += S.U;
			V_T += S.V;
			state = DEAD;

			// TODO: Quantize x and y location, place stokes values in matrix
		}
    }

    /*
     * rejection() - choose alpha and beta angles for scattering
     */
	void rejection(void) {
        // Placeholder locals for choosing rejection angle alpha
        double theta, phi, Io, Iith;
        do {
            theta = acos(2*rand_num() - 1); // TODO: why do we need to do 2*rand - 1? isn't this redundant?
            phi = 2*M_PI*rand_num();
            Io = s11[0]*S.I + s12[0]*(S.Q*cos(2*phi) + S.U*sin(2*phi));
            int i = floor(theta*nangles / M_PI);
            Iith = s11[i]*S.I + s12[i]*(S.Q*cos(2*phi) + S.U*sin(2*phi));
        } while (rand_num()*Io <= Iith);

		alpha = theta;
		beta = phi;
		return;
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
     * scatter(void) - 
     * Performs scattering based on the polarizing angle and scattering
     * particles in the media. For spheres, use Mie scattering theory.
     */
    void scatter(void) {

        /* Step 1: Update directional cosine for alpha and beta */

        double cos_alpha = cos(alpha);
        double sin_alpha = sin(alpha);
        double cos_beta = cos(beta);
        double sin_beta = sin(beta); 
        double old_ui = U.i;
        double old_uj = U.j;
        double old_uk = U.k;

        // Shortcut for when we are just about perpendicular to z axis
        if (1 - U.k <= 1e-12) {
            U.i = sin_alpha*cos_beta;
            U.j = sin_alpha*sin_beta;
            if (U.k > 0) {
                U.k = cos_alpha;
            } else {
                U.k = -cos_alpha;
            }
        }
        else {
            // Placeholders for the calculation
            double sin_z = sqrt(1-(U.k)*(U.k));

            U.i = sin_alpha*(old_ui*old_uk*cos_beta - old_uj*sin_beta) / (sin_z + old_ui*cos_alpha);
            U.j = sin_alpha*(old_uj*old_uk*cos_beta + old_ui*sin_beta) / (sin_z + old_uj*cos_alpha);
            U.k = -sin_alpha*cos_beta*sin_z + old_uk*cos_alpha;
        }

        /* Step 2: Rotate Stokes vector by beta */
        // Placeholders for step 4
        double si = S.I;
        double sq = S.Q;
        double su = S.U;
        double sv = S.V;
        S.rotate_stokes(beta);

        /* Step 3: Multiply Stokes by scattering matrix */

        int ithedeg = floor(alpha*nangles/M_PI);

		S.I= s11[ithedeg]*si+s12[ithedeg]*sq;
			
		S.Q= s12[ithedeg]*si+s11[ithedeg]*sq;

		S.U= s33[ithedeg]*su+s43[ithedeg]*sv;
			
		S.V= -s43[ithedeg]*su+s33[ithedeg]*sv;
	    	
        /* Step 4: Rotate Stokes again to get back into reference frame */
        /* TODO: Figure out why this is a necessary step. */

        double denominator = sqrt(1-cos_alpha*cos_alpha)*sqrt(1-U.k*U.k);
        double cosi;
        if (denominator <= 1e-12) { cosi = 0; }
        else {
            if ((beta > M_PI) && (beta < 2*M_PI)) {
                cosi = (U.k*cos_alpha - old_uk) / denominator;
            } else {
                cosi = -(U.k*cos_alpha - old_uk) / denominator;
                if (cosi > 1) { cosi = 1; }
                else if (cosi < -1) { cosi = -1;}
            }
        }
        // Rotate stokes vector clockwise by i
        S.rotate_stokes(2*M_PI - acos(cosi));
        S.normalize();
    }
};

#endif
