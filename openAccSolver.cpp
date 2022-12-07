#include "openAccSolver.h"
#include "mpi.h"
#include <vector>
#include <array>
#include <cstdlib>

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

OpenAccSolver::OpenAccSolver(double l[3], int dots) {
    Lx = l[0];
    Ly = l[1];
    Lz = l[2];

    N = dots;

    hx = Lx / dots;
    hy = Ly / dots;
    hz = Lz / dots;

    tau = 1.0 / 1000;
    alpha = M_PI * sqrt((4/(l[0]*l[0]) + 1/(l[1]*l[1]) + 1/(l[2]*l[2])));

    int periods[3] = {0,0,0};
    int plainRank;

    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &plainRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Dims_create(size, 3, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, true, &comminucator);
    MPI_Comm_rank(comminucator, &rank);
    MPI_Cart_coords(comminucator, rank, 3, procCoords);

    MPI_Cart_shift(comminucator, 0, 1, &xPrv, &xNxt);
    MPI_Cart_shift(comminucator, 1, 1, &yPrv, &yNxt);
    MPI_Cart_shift(comminucator, 2, 1, &zPrv, &zNxt);

    // zero and last coordinates get from other processes
    layerSize = (dots / dims[0] + (xPrv >= 0) + (xNxt >= 0) + (procCoords[0] == dims[0] - 1)) * 
            (dots / dims[1] + (yPrv >= 0) + (yNxt >= 0) + (procCoords[1] == dims[1] - 1)) *
            (dots / dims[2] + (zPrv >= 0) + (zNxt >= 0) +(procCoords[2] == dims[2] - 1));

    xStart = dots / dims[0] * procCoords[0] - (xPrv >= 0);
    xEnd = dots / dims[0] * (procCoords[0] + 1) + (xNxt >= 0) + (procCoords[0] == dims[0] - 1);
    xSize = xEnd - xStart;
    yStart = dots / dims[1] * procCoords[1] - (yPrv >= 0);
    yEnd = dots / dims[1] * (procCoords[1] + 1) + (yNxt >= 0) + (procCoords[1] == dims[1] - 1);
    ySize = yEnd - yStart;
    zStart = dots / dims[2] * procCoords[2] - (zPrv >= 0);
    zEnd = dots / dims[2] * (procCoords[2] + 1) + (zNxt >= 0) + (procCoords[2] == dims[2] - 1);
    zSize = zEnd - zStart;

    // if exchange with prev => our data starts from 1
    innerXStart = (xPrv >= 0);
    innerYStart = (yPrv >= 0);
    innerZStart = (zPrv >= 0);
    // if echange with next => our data ends at len size-1
    innerXEnd = xSize - (xNxt >= 0);
    innerYEnd = ySize - (yNxt >= 0);
    innerZEnd = zSize - (zNxt >= 0);

}

double OpenAccSolver::laplasian(double *u_n, int i, int j, int k) {
    double u_ijk = u_n[getIndex(i,j,k)];
    double dx = (u_n[getIndex(i-1,j,k)] - 2*u_ijk + u_n[getIndex(i+1,j,k)]) / (hx * hx);
    double dy = (u_n[getIndex(i,j-1,k)] - 2*u_ijk + u_n[getIndex(i,j+1,k)]) / (hy * hy);
    double dz = (u_n[getIndex(i,j,k-1)] - 2*u_ijk + u_n[getIndex(i,j,k+1)]) / (hz * hz);
    return dx + dy + dz;
}

void OpenAccSolver::fill_borders(double *u_t, int t) {
    // periodic x
    
    for (int j = innerYStart; j < innerYEnd; j++) {
        for (int k = innerZStart; k < innerZEnd; k++) {
            double dot[3] = {0,(j+yStart)*hy,(k+zStart)*hz};
            double u_here = u_analytical(dot, tau*t);
            if (procCoords[0] == 0) u_t[getIndex(0,j,k)] = u_here;
            if (procCoords[0] == dims[0] - 1) u_t[getIndex(innerXEnd-1,j,k)] = u_here;
        }
    }
    // 1r y
    for (int i = innerXStart; i < innerXEnd; i++) {
        for (int k = innerZStart; k < innerZEnd; k++) {
            if (procCoords[1] == 0) u_t[getIndex(i,0,k)] = 0;
            if (procCoords[1] == dims[1] - 1) u_t[getIndex(i,innerYEnd-1,k)] = 0;
        }
    }
    // 1r z
    for (int i = innerXStart; i < innerXEnd; i++) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            if (procCoords[2] == 0) u_t[getIndex(i,j,0)] = 0;
            if (procCoords[2] == dims[2] - 1) u_t[getIndex(i,j,innerZEnd-1)] = 0;
        }
    }
}

void OpenAccSolver::fill_interface(double *u) {
    int xLen = innerXEnd - innerXStart;
    int yLen = innerYEnd - innerYStart;
    int zLen = innerZEnd - innerZStart;
    int t = 0;
    double xPrevRecv[yLen*zLen], xNextRecv[yLen*zLen], yPrevRecv[xLen*zLen], yNextRecv[xLen*zLen], zPrevRecv[yLen*xLen], zNextRecv[yLen*xLen];
    double xPrevSend[yLen*zLen], xNextSend[yLen*zLen], yPrevSend[xLen*zLen], yNextSend[xLen*zLen], zPrevSend[yLen*xLen], zNextSend[yLen*xLen];
    MPI_Request reqSend[6] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    MPI_Request reqRecv[6] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    MPI_Status stat[6];
    // make send/recv requests
    if (xPrv >= 0) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                xPrevSend[t] = u[getIndex(1,j,k)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(xPrevSend, yLen*zLen, MPI_DOUBLE, xPrv, rank, comminucator, &reqSend[0]);
        MPI_Irecv(xPrevRecv, yLen*zLen, MPI_DOUBLE, xPrv, xPrv, comminucator, &reqRecv[0]);
    }
    if (xNxt >= 0) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                xNextSend[t] = u[getIndex(innerXEnd-1,j,k)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(xNextSend, yLen*zLen, MPI_DOUBLE, xNxt, rank, comminucator, &reqSend[1]);
        MPI_Irecv(xNextRecv, yLen*zLen, MPI_DOUBLE, xNxt, xNxt, comminucator, &reqRecv[1]);
    }
    if (yPrv >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                yPrevSend[t] = u[getIndex(i,1,k)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(yPrevSend, xLen*zLen, MPI_DOUBLE, yPrv, rank, comminucator, &reqSend[2]);
        MPI_Irecv(yPrevRecv, xLen*zLen, MPI_DOUBLE, yPrv, yPrv, comminucator, &reqRecv[2]);
    }
    if (yNxt >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                yNextSend[t] = u[getIndex(i,innerYEnd-1,k)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(yNextSend, xLen*zLen, MPI_DOUBLE, yNxt, rank, comminucator, &reqSend[3]);
        MPI_Irecv(yNextRecv, xLen*zLen, MPI_DOUBLE, yNxt, yNxt, comminucator, &reqRecv[3]);
    }
    if (zPrv >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int j = innerYStart; j < innerYEnd; j++) {
                zPrevSend[t] = u[getIndex(i,j,1)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(zPrevSend, yLen*xLen, MPI_DOUBLE, zPrv, rank, comminucator, &reqSend[4]);
        MPI_Irecv(zPrevRecv, yLen*xLen, MPI_DOUBLE, zPrv, zPrv, comminucator, &reqRecv[4]);
    }
    if (zNxt >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int j = innerYStart; j < innerYEnd; j++) {
                zNextSend[t] = u[getIndex(i,j,innerZEnd-1)];
                t++;
            }
        }
        t = 0;
        MPI_Isend(zNextSend, yLen*xLen, MPI_DOUBLE, zNxt, rank, comminucator, &reqSend[5]);
        MPI_Irecv(zNextRecv, yLen*xLen, MPI_DOUBLE, zNxt, zNxt, comminucator, &reqRecv[5]);
    }

    // fill u with recieved values
    MPI_Waitall(6, reqRecv, stat);
    if (xPrv >= 0) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                u[getIndex(0,j,k)] = xPrevRecv[t];
                t++;
            }
        }
        t = 0;
    }
    if (xNxt >= 0) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                u[getIndex(xSize-1,j,k)] = xNextRecv[t];
                t++;
            }
        }
        t = 0;
    }
    if (yPrv >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                u[getIndex(i,0,k)] = yPrevRecv[t];
                t++;
            }
        }
        t = 0;
    }
    if (yNxt >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                u[getIndex(i,ySize-1,k)] = yNextRecv[t];
                t++;
            }
        }
        t = 0;
    }
    if (zPrv >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int j = innerYStart; j < innerYEnd; j++) {
                u[getIndex(i,j,0)] = zPrevRecv[t];
                t++;
            }
        }
        t = 0;
    }
    if (zNxt >= 0) {
        for (int i = innerXStart; i < innerXEnd; i++) {
            for (int j = innerYStart; j < innerYEnd; j++) {
                u[getIndex(i,j,zSize-1)] = zNextRecv[t];
                t++;
            }
        }
        t = 0;
    }
    MPI_Waitall(6, reqSend, stat);

}

double OpenAccSolver::find_error(double *u_t, int t) {
    double eps, max_eps = 0;
    for (int i = innerXStart; i < innerXEnd; i++) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                double dot[3] = {hx*(i+xStart), hy*(j+yStart), hz*(k+zStart)};
                eps = fabs(u_t[getIndex(i,j,k)] - u_analytical(dot, tau*t));
                if (eps > max_eps) {
                    max_eps = eps;
                }
            }
        }
    }
    return max_eps;
}

double OpenAccSolver::solve(int steps) {
    double u[steps][layerSize];
    double error, maxError;

    for (int i = innerXStart; i < innerXEnd; i++) {
        for (int j = innerYStart; j < innerYEnd; j++) {
            for (int k = innerZStart; k < innerZEnd; k++) {
                double dot[3] = {hx*(i+xStart), hy*(j+yStart), hz*(k+zStart)};
                u[0][getIndex(i,j,k)] = phi(dot);
            }
        }
    }

    error = find_error(u[0], 0);
    MPI_Reduce(&error, &maxError, 1, MPI_DOUBLE, MPI_MAX, 0, comminucator);
    if (!rank) std::cout << maxError << std::endl;

    fill_interface(u[0]);
    fill_borders(u[1], 1);
    for (int i = 1; i < xSize-1; i++) {
        for (int j = 1; j < ySize-1; j++) {
            for (int k = 1; k < zSize-1; k++) {
                u[1][getIndex(i,j,k)] = u[0][getIndex(i,j,k)] + laplasian(u[0], i, j, k)*tau*tau/2;
            }
        }
    }

    error = find_error(u[1], 1);
    MPI_Reduce(&error, &maxError, 1, MPI_DOUBLE, MPI_MAX, 0, comminucator);
    if (!rank) std::cout << maxError << std::endl;

    for (int t = 2; t < steps; t++) {
        fill_interface(u[t-1]);
        fill_borders(u[t], t);
        for (int i = 1; i < xSize-1; i++) {
            for (int j = 1; j < ySize-1; j++) {
                for (int k = 1; k < zSize-1; k++) {
                    u[t][getIndex(i,j,k)] = 2 * u[t-1][getIndex(i,j,k)] - u[t-2][getIndex(i,j,k)] + tau * tau * laplasian(u[t-1], i, j, k);
                }
            }
        }

        error = find_error(u[t], t);
        MPI_Reduce(&error, &maxError, 1, MPI_DOUBLE, MPI_MAX, 0, comminucator);
        if (!rank) std::cout << maxError << std::endl;
        MPI_Barrier(comminucator);
    }

    return maxError;
}