#include <iostream>
#include <locale>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h>
#include <vector>
#include <random>

using namespace std;

#include "MPIUtilities.h"

//#define MPI_Matrix_multiply_vectors;

void mpiMatrixMultiply(int size, int rank);
vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int arows, int acols, int brows, int bcols, int procCount, int rank);
vector<int> mpiSort(vector<int> a, int size, int rank);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);

int MIN = 0;
int MAX = 100;
const double PI = 3.141592;

// 0 -> Pi/6
double F18(double x) {
    return x * cos(3 * x);
}

#pragma region MPI_Matrix_multiply_variables
const int arows = 500;
const int acols = 500;
const int brows = 500;
const int bcols = 500;
const int rowsPerProccess = 125; // rowsPerProccess = arows/size
static int a[arows][acols];
static int b[brows][bcols];
int c[arows][bcols];
#pragma endregion
#pragma region MPI_Sort_variables
const int sortNumsCount = 5000000;
#pragma endregion
#pragma region MPI_Math_variables

#pragma endregion

int main(int argc, char* argv[])
{
    double start, end; // Timer
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#pragma region MPI_Matrix_multiply
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
        cout << "Matrix mult "<< arows << "x" << acols << " * " << brows << "x" << bcols << " - " << abs(start - end) << " sec" << endl;
    }
#pragma endregion
#pragma region MPI_Matrx_multiply_vectors
#if defined MPI_Matrix_multiply_vectors
    vector<vector<int>> a(8, vector<int>(8));
    vector<vector<int>> b(8, vector<int>(8));
    vector<vector<int>> c(8, vector<int>(8));

    if (rank == 0)
    {
        a = vectorRandom(8, 8);
        b = vectorRandom(8, 8);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        start = MPI_Wtime();

    c = mpiMatrixMultiply(a, b, 8, 8, 8, 8, size, rank);
    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Matrix mult " << 8 << "x" << 8 << " * " << 8 << "x" << 8 << " - " << abs(start - end) << " sec" << endl;
    }
#endif // MPI_Matrix_multiply_vectors
#pragma endregion
#pragma region MPI_Sort
    static vector<int> t(sortNumsCount);
    
    if (rank == 0) {
        t = vectorRandom(sortNumsCount);
        start = MPI_Wtime();
    }

    auto c1 = mpiSort(t, size, rank);
    
    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Sort " << sortNumsCount << " - " << abs(start - end) << " sec" << endl;
    }
#pragma endregion
#pragma region MPI_Math

#pragma endregion

    MPI_Finalize();
}
vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int arows, int acols, int brows, int bcols, int procCount, int rank) {
    if (arows % procCount != 0)
        throw new exception("Количество строчек должно быть кратно количеству процессов");
    if (acols != brows)
        throw new exception("Нельзя перемножить");

    int sum = 0;
    int rowsPerProccess = arows / procCount;
    int resultSize = arows * bcols / procCount;
    
    vector<vector<int>> aa(rowsPerProccess, vector<int>(acols));
    vector<vector<int>> cc(rowsPerProccess, vector<int>(bcols));
    vector<vector<int>> c(arows, vector<int>(bcols));
    MPI_Scatter(&a.front(), rowsPerProccess * acols, MPI_INT, &aa.front(), rowsPerProccess * acols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b.front(), brows * bcols, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 1)
        cout << b[0][1] << endl;
    
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < rowsPerProccess; i++)
        for (int j = 0; j < bcols; j++) {
            for (int s = 0; s < acols; s++)
                sum += aa[i][s] * b[s][j];
            cc[i][j] = sum;
            sum = 0;
        }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&cc.front(), resultSize, MPI_INT, &c.front(), resultSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return c;
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
    auto aa = new int[rowsPerProccess][acols];
    auto cc = new int[rowsPerProccess][bcols];
    
    MPI_Barrier(MPI_COMM_WORLD);
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
vector<int> mpiSort(vector<int> a, int size, int rank)
{
    auto numsPerProccess = a.size() / size;
    static vector<int> c(a.size());
    static vector<int> aa(numsPerProccess);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Scatter(&a.front(), numsPerProccess, MPI_INT, &aa.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    sort(aa.begin(), aa.end());
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(&aa.front(), numsPerProccess, MPI_INT, &c.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        sort(c.begin(), c.end());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return c;
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