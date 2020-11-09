#include <iostream>
#include <cmath>


double estimationError (double * u, double * y, int N);
void realSolv(double * u, int N, double x0);
void sweepMethod(double * y, int N);


int main() {
    int N;
    double eps = 0.01, x0 = 0.5;

    std::cout << "Please, enter N:" << std::endl;
    std::cin >> N;

    while (N < 100000) {
        double * u;
        double * y;

        u = new double [N + 1];
        y = new double [N + 1];

        realSolv(u, N, x0);
        sweepMethod(y, N);

        for (int i = 0; i < N + 1; i++) {
            std::cout << "u[i] = " << u[i] << ", y[i] = " << y[i] << std::endl;
        }

        if (estimationError(u, y, N) < eps) {
            break;
        } else {
            N = 2 * N;
        }

        delete [] u;
        delete [] y;
    }
    return 0;
}


double estimationError(double * u, double * y, int N) {
    double res;

    for (int i = 0; i < N + 1; i++) {
        if (res < abs(u[2 * i] - y[i]) / 3) {
            res = abs(u[2 * i] - y[i]);
        }
    }

    return res;
}


void realSolv(double * u, int N, double x0) {
    double * x;

    double h = 1 / N;
    
    double k = x0 * x0 + 1;
    double q = x0;
    double f = exp(-x0); 

    double C = f / (q * ((sqrt(k * q) - 1) * exp(-sqrt(q / k)) - (sqrt(q / k) + 1) * exp(sqrt(q / k))));
    x = new double [N + 1];

    x[0] = 0.0;
    for (int i = 1; i < N + 1; i++) {
        x[i] = i * h;
    }

    for (int i = 0; i < N + 1; i++) {
        u[i] = C * (exp(sqrt(q / k) * x[i]) + exp(-sqrt(q / k) * x[i])) + f / q;
    }

    delete [] x;
}


void sweepMethod(double * y, int N) {
    double * k;
    double * q;
    double * fi;
    double * x;
    
    double * alf;
    double * bet;
    double * a;
    double * b;
    double * c;

    double h = 1.0 / N;

    k = new double [N + 1];
    q = new double [N + 1];
    fi = new double [N + 1];
    x = new double [N + 1];

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

    delete [] k;
    delete [] q;
    delete [] fi;
    delete [] x;

    delete [] alf;
    delete [] bet;
    delete [] a;
    delete [] b;
    delete [] c;
}