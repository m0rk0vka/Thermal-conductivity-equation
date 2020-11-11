#include <iostream>
#include <cmath>


double estimationError (double * u, double * y, int N);
void realSolv(double * u, int N, double x0);
void sweepMethod(double * y, int N);


int main() {
    double * u;
    double * y;
    double * y2;

    int N;
    double eps = 0.01, x0 = 0.5;

    std::cout << "Please, enter N: ";
    std::cin >> N;

    int j = 0;

    while (N < 100000) {
        u = new double [N + 1];
        y = new double [N + 1];

        realSolv(u, N, x0);
        sweepMethod(y, N);

        if (j == 1) {
            if (estimationError(y2, y, N) < eps) {
                for (int i = 0; i < N + 1; i++) {
                    std::cout << "x[" << i << "] = " << i * 1.0 / N << " , u[" << i << "] = " << u[i] << " , y[" <<
                        i << "] = " << y[i] << " , delta = " << abs(u[i] - y[i]) << std::endl;
                }

                break;
            } else {
                delete [] y2;

                y2 = new double [N + 1];

                for (int i = 0; i < N + 1; i++) {
                    y2[i] = y[i];
                }

                N = 2 * N;
            }
        } else {
            y2 = new double [N + 1];

            for (int i = 0; i < N + 1; i++) {
                y2[i] = y[i];
            }

            N = 2 * N;
        }

        delete [] u;
        delete [] y;

        j = 1;
    }
    return 0;
}


double estimationError(double * y2, double * y, int N) {
    double res;

    N = N / 2;

    for (int i = 0; i < N + 1; i++) {
        if (res < abs(y[2 * i] - y2[i]) / 3) {
            res = abs(y[2 * i] - y2[i]) / 3;
        }
    }

    return res;
}


void realSolv(double * u, int N, double x0) {
    double * x;

    double h = 1.0 / N;
    
    double k = x0 * x0 + 1;
    double q = x0;
    double f = exp(-x0); 

    double C = f / (q * ((sqrt(k * q) - 1) * exp(-sqrt(q / k)) - (sqrt(k * q) + 1) * exp(sqrt(q / k))));
    x = new double [N + 1];

    for (int i = 0; i < N + 1; i++) {
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
        x[i] = i * h;

        /*k[i] = x[i] * x[i] + 1;
        q[i] = x[i];
        fi[i] = exp(-x[i]);*/

        k[i] = 1.25;
        q[i] = 0.5;
        fi[i] = exp(-0.5);
    }

    for (int i = 1; i < N; i++) {
        a[i] = 1 / (2 * h * h) * (k[i - 1] + k[i]);
        b[i] = 1 / (2 * h * h) * (k[i] + k[i + 1]);
        c[i] = q[i] + 1 / (2 * h * h) * (k[i - 1] + 2 * k[i] + k[i + 1]);
    }

    alf[1] = (k[0] + k[1]) / (k[0] + k[1] + h * h * q[0]);
    bet[1] = h * h * fi[0] / (k[0] + k[1] + h * h * q[0]);

    for (int i = 1; i < N; i++) {
        alf[i + 1] = b[i] / (c[i] - a[i] * alf[i]);
        bet[i + 1] = (fi[i] + a[i] * bet[i]) / (c[i] - a[i] * alf[i]); 
    }

    double mu1 = 1 + 0.5 * h * q[N] + (k[N] + k[N - 1]) / (2 * h);
    double hi2 = (k[N] + k[N - 1]) / (2 * h * mu1);

    y[N] = (0.5 * h * fi[N] / mu1 + hi2 * bet[N]) / (1 - alf[N] * hi2);

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