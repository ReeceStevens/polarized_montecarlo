/*
 * Vector.h
 * Physical <i,j,k> vector class for Monte Carlo simulation
 * along with some useful utility member functions
 */

#include <stdio.h>
#include <math.h>

struct Vector {
public:
    double i;
    double j;
    double k;
    Vector(void) { i = 0; j = 0; k = 0; };

    Vector(double i, double j, double k) : i(i), j(j), k(k){};

    Vector(const Vector& that) {
        this->i = that.i;
        this->j = that.j;
        this->k = that.k;
    };

    /* Mathematical Operations */

    /*
     * Vector::cross_prod(a,b) 
     * Return the cross product of vectors
     * a and b in new vector c.
     */
    static Vector cross_prod(Vector& a, Vector& b) {
        double new_i = a.j * b.k - a.k * b.j;
        double new_j = -(a.i * b.k - a.k * b.i);
        double new_k = a.i * b.j - a.j * b.i;
        Vector result(new_i, new_j, new_k);
        return result;
    };

    /*
     * Vector::dot_prod(a,b) 
     * Return the dot product of vectors
     * a and b.
     */
    static double dot_prod(Vector& a, Vector& b) {
        return a.i*b.i + a.j*b.j + a.k*b.k;
    };

    static Vector scalar_mult(double x, Vector &a) { 
        Vector result(a.i*x, a.j*x, a.k*x);
        return result;
    }

    // Three-element sum
    static Vector sum(Vector& a, Vector& b, Vector & c) {
        Vector result(a.i+b.i+c.i, a.j+b.j+c.j, a.k+b.k+c.k);
        return result;
    }

    void operator=(Vector& a) {
        i = a.i;
        j = a.j;
        k = a.k;
    }
};
