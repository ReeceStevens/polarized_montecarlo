#include <stdio.h>
#include <math.h>

struct Matrix {
public:
	int length;
	int width;

	Matrix(int rows, int columns){
		// Create the space.
		double** data = new double* [columns];
		for (int i = 0; i < columns; i += 1) {
			data[i] = new double [rows];
		}

		
	}
}
