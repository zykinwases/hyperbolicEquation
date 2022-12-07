#include <iostream>
#include <fstream>
#include <cstdlib>
#include "ompMpiSolver.h"
#include "mpiSolver.h"
#include "solver.h"
#include "mpi.h"

int main(int argc, char *argv[])
{
    int dots, solverMode;
    double L;
    std::ofstream output;

    // parse arguments
    // arguments -- solver, dots, L (assume Lx=Ly=Lz) (1 == 1, 2 == Pi), fname
    if (argc > 3) {
        solverMode = std::atoi(argv[1]);
        dots = std::atoi(argv[2]);
        int c = std::atoi(argv[3]);
        if (c == 2) L = M_PI;
        else L = 1;
        if (dots <= 0) dots = 2;
        if (solverMode <= 0) solverMode = 1;
        output.open(argv[4], std::ios_base::app);
    } else {
        solverMode = 1;
        dots = 2;
        L = 1;
    }

    double l[3] = {L,L,L};
    double error;
    double begin, time, maxTime;
    int size,rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    switch (solverMode) {
    case 1: {
        Solver solver(l, dots);
        begin = MPI_Wtime();
        error = solver.solve(20);
        time = MPI_Wtime() - begin;
        if (!rank) output << dots << " " << size << " " << error << " " << time << std::endl;
        break;
    }
    case 2: {
        MpiSolver solver(l, dots);
        begin = MPI_Wtime();
        error = solver.solve(20);
        time = MPI_Wtime() - begin;
        MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (!rank) output << dots << " " << size << " " << error << " " << maxTime << std::endl;
        break;
    }
    case 3: {
        OmpMpiSolver solver(l, dots);
        begin = MPI_Wtime();
        error = solver.solve(20);
        time = MPI_Wtime() - begin;
        MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (!rank) output << dots << " " << size << " " << error << " " << maxTime << std::endl;
    }
    }

    MPI_Finalize();

    output.close();

    return 0;
}
