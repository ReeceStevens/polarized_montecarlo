#include <stdio.h>
#include <assert.h>

class Complex{
private:
	double r;
	double i;

public:
	Complex(double r, double i) {
		this->r = r;
		this->i = i;		
	}

	Complex(const Complex& that) {
		this->r = that.r;
		this->i = that.i;
	}

	double Re() {
		return r;
	}

	double Im() {
		return i;
	}

	static Complex& conj(Complex a) {
		//Complex* retval = new Complex(this->r,-1*this->i);
		//Complex& ret = *retval;
		Complex& retval = *(new Complex(a.r,-1*a.i));
		return retval;
	}

	Complex& operator+(const Complex& that) {
		double real = this->r + that.r;
		double imag = this->i + that.i;
		Complex* result = new Complex(real,imag);
		Complex& ret = *result;
		return ret;	
	}

	bool operator==(const Complex& that) {
		return ((this->r == that.r) && (this->i == that.i));
	}

	Complex& operator*(const Complex& that) {
		double real = this->r*that.r - this->i*that.i;
		double imag = this->r*that.i + this->i*that.r;
		Complex& retval = *(new Complex(real,imag));
		return retval;		
	}

	Complex& operator/(const double k) {
		assert(k != 0);
		Complex& retval = *(new Complex(r/k,i/k));
		return retval; 
	}

	Complex& operator/(const Complex& that) {
		assert(((that.r != 0) || (that.i != 0)));
		return ((*this)*(conj(that)))/(that.r*that.r + that.i*that.i);
	}

};
