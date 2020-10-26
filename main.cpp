#include <iostream>

/* sweepMethod from wiki:
	a[0] = 0; b[n-1] = 0;	
*/
void sweepMethod(int n, double *a, double* b, double *c, double *f, double *x);


int main() {
	return 0;
}


void sweepMethod(int n, double *a, double* b, double *c, double *f, double *x) {
	
	double m;
	
	for (int i = 1; i < n; i++) {
		m = a[i] / c[i - 1];
		c[i] = c[i] - m * b[i - 1];
		f[i] = f[i] - m * f[i - 1];
	}

	x[n - 1] = f[n - 1] / c[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
	}
}