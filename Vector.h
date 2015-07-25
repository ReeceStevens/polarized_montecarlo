#include <stdio.h>
#include <math.h>

struct Vector {
public:
    double i;
    double j;
    double k;
    Vector(void) { i = 0; j = 0; k = 0; };

    Vector(double i, double j, double k) : i(i), j(j), k(k){};

    Vector(Vector& that) {
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
    Vector cross_prod(Vector& a, Vector& b) {
        double new_i = a.j * b.k + a.k * b.j;
        double new_j = a.i * b.k + a.k * b.i;
        double new_k = a.i * b.j + a.j * b.i;
        Vector result(new_i, new_j, new_k);   
        return result;
    };

    /*
     * Vector::dot_prod(a,b) 
     * Return the dot product of vectors
     * a and b.
     */
    double dot_prod(Vector& a, Vector& b) {
        return a.i*b.i + a.j*b.j + a.k*b.k;
    };
};
