#include <iostream>

/* sweepMethod from wiki:
    a[0] = 0; b[n-1] = 0;
*/
void sweepMethod(int N = 1000);


int main() {
    return 0;
}

//error estimate for the solution of a boundary value problem by Runge's rule
double estimationError(double *y, double *u, int N, int m) {
    
}

void sweepMethod(int N = 1000) {
    double * k, q, fi, x, y, u;
    double * alf, bet, a, b, c;
    double h = 1.0 / N;

    k = new double [N + 1];
    q = new double [N + 1];
    fi = new double [N + 1];
    x = new double [N + 1];
    y = new double [N + 1];
    u = new double [N + 1];

    alf = new double [N + 2];
    bet = new double [N + 2];
    a = new double [N + 2];
    b = new double [N + 2];
    c = new double [N + 2];

    for (int i = 0; i < N + 1; i++) {
        k[i] = q[i] = fi[i] = x[i] = 0.0;
    }

    for (int i = 1; i < N + 1; i++) {
        x[i] = i * h;

        k[i] = x[i] * x[i] + 1;
        q[i] = x[i];
        fi[i] = exp(-x[i]);
    }

    for (int i = 1; i < N; i++) {
        a[i] = 1 / (2 * h * h) * (k[i - 1] + k[i]);
        b[i] = 1 / (2 * h * h) * (k[i] + k[i + 1]);
        c[i] = q[i] + 1 / (2 * h * h) * (k[i - 1] + 2 * k[i] + k[i + 1]);
    }

    alf[1] = 0.0; // b[0] / c[0];
    bet[1] = 0.0; // fi[0] / c[0];

    for (int i = 1; i < N; i++) {
        alf[i + 1] = b[i] / (c[i] - a[i] * alf[i]);
        bet[i + 1] = (fi[i] + alf[i] * bet[i]) / (c[i] - a[i] * alf[i]); 
    }

    y[N] = -((1 + 0.5 * h * q[N]) * 1 - 0.5 * h * fi[N]);

    for (int i = N - 1; i >= 0; i--) {
        y[i] = alf[i + 1] * y[i + 1] + bet[i + 1];
    }
}