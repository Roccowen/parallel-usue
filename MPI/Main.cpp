#include <iostream>
#include <locale>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h>

using namespace std;

#include "MPITest.h"
#include "MPIUtilities.h"

void mpiMatrixMultiply(int size, int rank);

int MIN = 0;
int MAX = 10;
const double PI = 3.141592;

const int arows = 500;
const int acols = 500;
const int brows = 500;
const int bcols = 500;
const int rowsPerProccess = 125; // arows/size

static int a[arows][brows];
static int b[brows][bcols];
int c[arows][bcols];

int main(int argc, char* argv[])
{
    double start, end;
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        for (size_t i = 0; i < arows; i++)
            for (size_t j = 0; j < acols; j++)
                a[i][j] = randomInt();
        for (size_t i = 0; i < brows; i++)
            for (size_t j = 0; j < bcols; j++)
                b[i][j] = randomInt();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        start = MPI_Wtime();
    mpiMatrixMultiply(size, rank);
    if (rank == 0) {
        end = MPI_Wtime();
        cout << abs(start - end) << endl;
    }


    MPI_Finalize();
}

void mpiMatrixMultiply(int size, int rank) {
    if (arows % size != 0)
        throw new exception("Количество строчек должно быть кратно количеству процессов");
    if (acols != brows)
        throw new exception("Нельзя перемножить");

    int sum = 0;
    int matrixSize = acols;    
    int partSize = arows * rowsPerProccess;
    int resultSize = arows * bcols / size;
    MPI_Barrier(MPI_COMM_WORLD);

    auto aa = new int[rowsPerProccess][acols];
    auto cc = new int[rowsPerProccess][bcols];

    MPI_Scatter(a, partSize, MPI_INT, aa, partSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, brows * bcols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < rowsPerProccess; i++)
        for (int j = 0; j < bcols; j++) {
            for (int s = 0; s < matrixSize; s++)
                sum += aa[i][s] * b[s][j];
            cc[i][j] = sum;
            sum = 0;
        }
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(cc, resultSize, MPI_INT, c, resultSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}

//void mpiSort(vector<int>& arr)
//{
//    int length = arr.size();
//    int* arrayPtr = &arr[0];
//    int temp,
//        item;
//    for (int counter = 1; counter < length; counter++)
//    {
//        temp = arrayPtr[counter];
//        item = counter - 1;
//        while (item >= 0 && arrayPtr[item] > temp)
//        {
//            arrayPtr[item + 1] = arrayPtr[item];
//            arrayPtr[item] = temp;
//            item--;
//        }
//    }
//}
//double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n) {
//    double step = (b - a) / n;
//    double result = 0;
//    for (double x = a + step; x < b; x += step)
//        result += F(x);
//    return (result + F(a) / 2 + F(b) / 2) * step;
//}
//double mpiRectangleMethod(double (*F)(double), double a, double b, int n) {
//    double step = (b - a) / n;
//    double result = 0;
//    for (double x = a + step; x <= b; x += step)
//        result += F(x);
//    return result * step;
//}