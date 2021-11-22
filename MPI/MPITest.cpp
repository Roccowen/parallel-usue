#include <vector>
#include <random>
#include <mpi.h>
#include <vector>

using namespace std;

#include "MPIUtilities.h"
#include "MPITest.h"
#include <iostream>

int** mpiMatrixMultiply(int** a, int** b, int aRows, int aCols, int bRows, int bCols, int size)
{
    if (size != aRows)
        throw new exception("���������� ��������� ������� � ������� ������� ������ ���������!");

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int** c = new int* [aRows];
    for (int i = 0; i < aRows; ++i)
        c[i] = new int[bCols];

    auto aa = new int[aCols];
    auto cc = new int[bCols];
    int  sum = 0;

    MPI_Scatter(a, aCols, MPI_INT, aa, aCols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, bRows * bCols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < aCols; i++)
    {
        cout << "5 - " << rank << endl;
        for (int j = 0; j < bRows; j++) {
            cout << "6 - " << rank << endl;
            sum += aa[j] * b[j][i];
        }
        cout << "7 - " << rank << endl;
        cc[i] = sum;
        sum = 0;
    }
    
    MPI_Gather(cc, aRows * bCols / size, MPI_INT, c, aRows * bCols / size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "8" << endl;
    return c;
}
void mpiMatrixMultiplyTest(int size) {

    const int arows = 8;
    const int acols = 8;
    const int brows = 8;
    const int bcols = 8;
    
    if (arows % size != 0)
        throw new exception("����� �����-��, �����������");

    int sum = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a[arows][acols];
    int b[brows][bcols];

    int partSize = arows * acols / size;
    int partRowsCount = arows / size;

    if (rank == 0)
    {
        cout << "A" << endl;
        for (size_t i = 0; i < arows; i++) {
            for (size_t j = 0; j < acols; j++) {           
                a[i][j] = randomInt();
                cout << a[i][j] << " ";
            }
            cout << endl;
        }

        cout << "B" << endl;
        for (size_t i = 0; i < brows; i++) {
            for (size_t j = 0; j < bcols; j++) {
                b[i][j] = randomInt();
                cout << b[i][j] << " ";
            }
            cout << endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int c[arows][bcols];
    auto aa = new int[partRowsCount][acols];
    auto cc = new int[partSize];
   
    MPI_Scatter(a, partSize, MPI_INT, aa, partSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, brows * bcols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
    {
        cout << endl;
        for (size_t i = 0; i < partRowsCount; i++) {
            for (size_t j = 0; j < acols; j++) {
                cout << aa[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    for (int i = 0; i < partRowsCount; i++)
    {
        for (int j = 0; j < acols; j++)
        {
            for (int s = 0; s < brows; s++)
                //sum += aa[j + (acols * i)] * b[j][i];
            cc[i + (acols * i)] = sum;
            sum = 0;
        }
    }


    MPI_Gather(cc, partSize, MPI_INT, c, partSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << "C" << endl;
        for (size_t i = 0; i < arows; i++) {
            for (size_t j = 0; j < bcols; j++) {
                cout << c[i][j] << " ";
            }
            cout << endl;
        }
    }

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