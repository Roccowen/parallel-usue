﻿#include <iostream>
#include <locale>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h>

using namespace std;

#include "MPITest.h"

void MPITest();


int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    MPI_Finalize(); 
}



void MPITest() {

}

