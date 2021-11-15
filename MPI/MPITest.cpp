#include <vector>
#include <random>
#include <mpi.h>
#include <vector>

using namespace std;

#include "MPITest.h"

int** mpiMatrixMultiply(int** a, int** b, int aRows, int aCols, int bRows, int bCols, int size)
{
    auto aLocal = new int[aRows * aCols / size];
    auto cLocal = new int[aRows * bCols / size];
    int temp;
    
    int** ñ = new int* [aRows];
    for (int i = 0; i < aRows; ++i)
        ñ[i] = new int[bCols];
   
    MPI_Scatter(a, aRows * aCols / size, MPI_INT, aLocal, aRows * aCols / size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, bRows * bCols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < bRows; i++)
    {
        for (int j = 0; j < bCols; j++)
            temp += aLocal[j] * b[j][i];
        cLocal[i] = temp;
        temp = 0;
    }

    MPI_Gather(cLocal, aRows * bCols / size, MPI_INT, ñ, aRows * bCols / size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return ñ;
}
void mpiSort(vector<int>& arr)
{
    int length = arr.size();
    int* arrayPtr = &arr[0];
    int temp,
        item;
    for (int counter = 1; counter < length; counter++)
    {
        temp = arrayPtr[counter];
        item = counter - 1;
        while (item >= 0 && arrayPtr[item] > temp)
        {
            arrayPtr[item + 1] = arrayPtr[item];
            arrayPtr[item] = temp;
            item--;
        }
    }
}
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double result = 0;
    for (double x = a + step; x < b; x += step)
        result += F(x);
    return (result + F(a) / 2 + F(b) / 2) * step;
}
double mpiRectangleMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double result = 0;
    for (double x = a + step; x <= b; x += step)
        result += F(x);
    return result * step;
}