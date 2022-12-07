#include "solver.h"
#include <vector>
#include <array>
#include <cstdlib>

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

Solver::Solver(double l[3], int dots) {
    Lx = l[0];
    Ly = l[1];
    Lz = l[2];

    N = dots;
    layerSize = (dots+1)*(dots+1)*(dots+1);
    hx = Lx / dots;
    hy = Ly / dots;
    hz = Lz / dots;

    tau = 1.0 / 1000;
    alpha = M_PI * sqrt((4/(l[0]*l[0]) + 1/(l[1]*l[1]) + 1/(l[2]*l[2])));
}

double Solver::laplasian(double *u_n, int i, int j, int k) {
    double u_ijk = u_n[getIndex(i,j,k)];
    double dx = (u_n[getIndex(i-1,j,k)] - 2*u_ijk + u_n[getIndex(i+1,j,k)]) / (hx * hx);
    double dy = (u_n[getIndex(i,j-1,k)] - 2*u_ijk + u_n[getIndex(i,j+1,k)]) / (hy * hy);
    double dz = (u_n[getIndex(i,j,k-1)] - 2*u_ijk + u_n[getIndex(i,j,k+1)]) / (hz * hz);
    return dx + dy + dz;
}

void Solver::fill_borders(double *u_t, int t) {
    // periodic x
    for (int j = 0; j <= N; j++) {
        for (int k = 0; k <= N; k++) {
            double dot[3] = {0,j*hy,k*hz};
            double u_here = u_analytical(dot, tau*t);
            u_t[getIndex(0,j,k)] = u_here;
            u_t[getIndex(N,j,k)] = u_here;
        }
    }
    // 1r y
    for (int i = 0; i <= N; i++) {
        for (int k = 0; k <= N; k++) {
            u_t[getIndex(i,0,k)] = 0;
            u_t[getIndex(i,N,k)] = 0;
        }
    }
    // 1r z
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            u_t[getIndex(i,j,0)] = 0;
            u_t[getIndex(i,j,N)] = 0;
        }
    }
}

double Solver::find_error(double *u_t, int t) {
    double eps, max_eps = 0;
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {
                double dot[3] = {hx*i, hy*j, hz*k};
                eps = fabs(u_t[getIndex(i,j,k)] - u_analytical(dot, tau*t));
                if (eps > max_eps) {
                    max_eps = eps;
                }
            }
        }
    }
    return max_eps;
}

double Solver::solve(int steps) {
    double u[steps][layerSize];

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {
                double dot[3] = {hx*i, hy*j, hz*k};
                u[0][getIndex(i,j,k)] = phi(dot);
            }
        }
    }

    std::cout << find_error(u[0], 0) << std::endl;

    fill_borders(u[1], 1);
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            for (int k = 1; k < N; k++) {
                u[1][getIndex(i,j,k)] = u[0][getIndex(i,j,k)] + laplasian(u[0], i, j, k)*tau*tau/2;
            }
        }
    }
    
    std::cout << find_error(u[1], 1) << std::endl;

    double error;
    for (int t = 2; t < steps; t++) {
        fill_borders(u[t], t);
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < N; j++) {
                for (int k = 1; k < N; k++) {
                    u[t][getIndex(i,j,k)] = 2 * u[t-1][getIndex(i,j,k)] - u[t-2][getIndex(i,j,k)] + tau * tau * laplasian(u[t-1], i, j, k);
                }
            }
        }

        error = find_error(u[t], t);
        std::cout << error << std::endl;
    }

    return error;
}
