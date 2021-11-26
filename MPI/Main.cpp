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

void mpiMatrixMultiply(int size, int rank);
vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int size, int rank);
vector<int> mpiSort(vector<int> a, int size, int rank);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);


int MIN = 0;
int MAX = 10;
const double PI = 3.141592;

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
#pragma region MPI_Sort_varizbles

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
    static vector<vector<int>> a;
    static vector<vector<int>> b;
    static vector<vector<int>> c;
    if (rank == 0)
    {
        a = vectorRandom(500, 500);
        b = vectorRandom(500, 500);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        start = MPI_Wtime();
    c = mpiMatrixMultiply(a, b, size, rank);
    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Matrix mult "<< arows << "x" << acols << " * " << brows << "x" << bcols << " - " << abs(start - end) << " sec" << endl;
    }
#pragma endregion
#pragma region MPI_Sort
    static vector<int> t(500);
    
    if (rank == 0)
    {
        t = vectorRandom(500);
    }
    auto c1 = mpiSort(t, size, rank);
    if (rank == 0) {
        end = MPI_Wtime();
        cout << "Sort " << 500 << " - " << abs(start - end) << " sec" << endl;
        printVector(c1);
    }
#pragma endregion
#pragma region MPI_Math

#pragma endregion

    MPI_Finalize();
}
vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int procCount, int rank) {
    if (a.size() == 0 || b.size() == 0)
        throw new exception("a и b должны иметь элементы");
    int rows = a.size();
    int size = a[0].size();
    int cols = b[0].size();
    if (rows % size != 0)
        throw new exception("Количество строчек должно быть кратно количеству процессов");
    if (a[0].size() != b.size())
        throw new exception("Нельзя перемножить");
    int sum = 0;
    int rowsPerProccess = rows / procCount;
    int resultSize = rows * cols / procCount;
    vector<vector<int>> aa(rowsPerProccess, vector<int>(size));
    vector<vector<int>> cc(rowsPerProccess, vector<int>(cols));
    vector<vector<int>> c(rows, vector<int>(cols));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(a.data(), rowsPerProccess, MPI_INT, aa.data(), rowsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b.data(), size * cols, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < rowsPerProccess; i++)
        for (int j = 0; j < cols; j++) {
            for (int s = 0; s < size; s++)
                sum += aa[i][s] * b[s][j];
            cc[i][j] = sum;
            sum = 0;
        }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(cc.data(), resultSize, MPI_INT, c.data(), resultSize, MPI_INT, 0, MPI_COMM_WORLD);
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
    
    MPI_Scatter(a.data(), numsPerProccess, MPI_INT, aa.data(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    sort(aa.begin(), aa.end());
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(aa.data(), numsPerProccess, MPI_INT, c.data(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
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