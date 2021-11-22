#include <iostream>
#include <locale>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h>

using namespace std;

#include "MPITest.h"
#include "MPIUtilities.h"

void MPITest(int argc, char* argv[]);

const double PI = 3.141592;
int MIN = 0;
int MAX = 10;

int main(int argc, char* argv[])
{

    MPITest(argc, argv);
}

void MPITest(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    mpiMatrixMultiplyTest(size);

    //if (rank == 0) {
    //    a = arrayRandom(4, 4);
    //    b = arrayRandom(4, 4);
    //    
    //    printArray(a, 4, 4);
    //    printArray(b, 4, 4);
    //}
    //
    //MPI_Barrier(MPI_COMM_WORLD);
    //
    //c = mpiMatrixMultiply(a, b, 4, 4, 4, 4, size);
    //
    //MPI_Barrier(MPI_COMM_WORLD);
    //
    //if (rank == 0) {
    //    printArray(arrayRandom(4, 4), 4, 4);
    //    auto a = arrayRandom(4, 4);
    //    auto b = arrayRandom(4, 4);
    //    printArray(c, 4, 1);
    //}


    MPI_Finalize();
}