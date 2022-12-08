#ifndef SOLVER
#define SOLVER

#include <vector>
#include <cmath>

class Solver {
private:
    int N, layerSize;
    double Lx, Ly, Lz;
    double hx, hy, hz;
    double tau, alpha;

    int getIndex(int x, int y, int z) {
        return (x*(N+1) + y)*(N+1) + z;
    };
    void fill_borders(double *u_t, int t);

    double u_analytical(double d[3], double t) {
        return (sin(2 * M_PI * d[0] / Lx) * sin(M_PI * d[1] / Ly) * sin(M_PI * d[2] / Lz) * cos(this->alpha * t + M_PI));
    };
    double phi(double d[3]) {
        // u_analytical(d, 0);
        return (-1 * sin(2 * M_PI * d[0] / Lx) * sin(M_PI * d[1] / Ly) * sin(M_PI * d[2] / Lz));
    };
    double laplasian(double *u_n, int i, int j, int k);

public:
    Solver(double l[3], int N);
    double solve(int i);
    double find_error(double *u_t, int t);
    ~Solver() {};
};

#endif