#ifndef OMPMPISOLVER
#define OMPMPISOLVER

#include <vector>
#include <cmath>
#include "mpi.h"

class OmpMpiSolver {
private:
    int N, layerSize;
    double Lx, Ly, Lz;
    double hx, hy, hz;
    double tau, alpha;
    int initialized = 0;

    MPI_Comm comminucator;
    int rank, size;
    int procCoords[3];
    int xPrv, xNxt, yPrv, yNxt, zPrv, zNxt;
    int xStart, xEnd, yStart, yEnd, zStart, zEnd;
    int xSize, ySize, zSize;
    int innerXStart, innerXEnd, innerYStart, innerYEnd, innerZStart, innerZEnd;
    int dims[3] = {0,0,0};

    int getIndex(int x, int y, int z) {
        return x + y*(xSize) + z*(xSize)*(ySize);
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
    void get_interface(double *u, double *xPrev, double *xNext, double *yPrev, double *yNext, double *zPrev, double *zNext);
    void fill_interface(double *u);

public:
    OmpMpiSolver(double l[3], int N);
    double solve(int i);
    double find_error(double *u_t, int t);
    ~OmpMpiSolver() {if (!initialized) MPI_Finalize();};
};

#endif