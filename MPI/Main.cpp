#include <iomanip>
#include <iostream>
#include <locale>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h>
#include <vector>
#include <random>
#include <chrono>

using namespace std;

#include "MPIUtilities.h"

#define MAX_SIZE 32764‬

#define MPI_Matrix_multiply
#define MPI_Sort
#define MPI_Rectangle_Method
#define MPI_Trapeziodal_Method

#ifdef MPI_Matrix_multiply
#define arows 500
#define acols 500
#define brows 500
#define bcols 500
// rowsPerProccess = arows/size.
#define rowsPerProccess 125  
static int a[arows][acols];
static int b[brows][bcols];
static int c[arows][bcols];
#endif
#ifdef MPI_Sort
#define sortNumsCount 5000000
#endif
#ifdef MPI_Rectangle_Method
#define stepsCount 1000
#endif
#ifdef MPI_Trapeziodal_Method
#define stepsCount 1000
#endif

void mpiMatrixMultiply(int a[arows][acols], int b[brows][bcols], int c[arows][bcols], int size, int rank);
vector<int> mpiSort(vector<int> a, int size, int rank);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n, int size, int rank);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n, int size, int rank);

int MIN = 0;
int MAX = 100;
const double PI = 3.141592;
// 0 -> Pi/6
double F18(double x) {
    return x * cos(3 * x);
}
int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "Russian");
    cout << std::setprecision(12);
    double start, end; // timer
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef MPI_Matrix_multiply
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
    if (rank == 0)
        start = MPI_Wtime();
    mpiMatrixMultiply(a, b, c, size, rank);
    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Matrix mult " << arows << "x" << acols << " * " << brows << "x" << bcols << " - " << abs(start - end) << " sec" << endl;
    }
#endif
#ifdef MPI_Sort
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
#endif
#ifdef MPI_Rectangle_Method
    if (rank == 0) {
        start = MPI_Wtime();
    }
    auto recResult = mpiRectangleMethod(F18, 0, PI / 6, stepsCount, size, rank);

    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Rectangle method " << stepsCount << " - " << abs(start - end) << " sec" << endl;
        cout << recResult << endl;
    }
#endif
#ifdef MPI_Trapeziodal_Method
    if (rank == 0) {
        start = MPI_Wtime();
    }
    auto trResult = mpiTrapeziodalMethod(F18, 0, PI / 6, stepsCount, size, rank);

    if (rank == 0)
    {
        end = MPI_Wtime();
        cout << "Trapeziodal method " << stepsCount << " - " << abs(start - end) << " sec" << endl;
        cout << trResult << endl;
    }
#endif
    MPI_Finalize();
}
void mpiMatrixMultiply(int a[arows][acols], int b[brows][bcols], int c[arows][bcols], int size, int rank) {
    MPI_Barrier(MPI_COMM_WORLD);

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
    MPI_Barrier(MPI_COMM_WORLD);
    auto numsPerProccess = a.size() / size;
    static vector<int> c(a.size());
    static vector<int> aa(numsPerProccess);
    
    MPI_Scatter(&a.front(), numsPerProccess, MPI_INT, &aa.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    sort(aa.begin(), aa.end());
    
    MPI_Gather(&aa.front(), numsPerProccess, MPI_INT, &c.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        sort(c.begin(), c.end());
    
    MPI_Barrier(MPI_COMM_WORLD);
    return c;
}
double mpiRectangleMethod(double (*F)(double), double a, double b, int n, int size, int rank) 
{ 
    MPI_Barrier(MPI_COMM_WORLD);
    if (n % size != 0)
        throw new exception("Количество шагов должно быть чётно количеству процессов");

    double step = (b - a) / n;
    double result = 0;
    double localResult = 0;
    for (double x = a + step * (rank + 1); x <= b; x += step * size)
        localResult += F(x);
    MPI_Reduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    return result * step;
}
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n, int size, int rank) 
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (n % size != 0)
        throw new exception("Количество шагов должно быть чётно количеству процессов");

    double step = (b - a) / n;
    double result = 0;
    double localResult = 0;
    
    for (double x = a + step * (rank + 1); x < b; x += step * size)
        localResult += F(x);
    
    MPI_Reduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return (result + F(a) / 2 + F(b) / 2) * step;
}