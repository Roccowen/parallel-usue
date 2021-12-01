#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <conio.h>
#include <locale.h>
#include <stdlib.h>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <time.h>
#include <functional>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <assert.h>


double f(double x)
{
    return (x + 1) * log(x);
}
void ParallelMatrixV()//500х500 * 500
{
    int rank, size;
    int i, j, n = 500;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n_partial = n / size;
    double* a_partial = new double[n_partial * n];//блоки строк исходной матрицы на каждом процессе
    double* x = new double[n]; //исходный вектор
    double* y_partial = new double[n_partial];//блоки результирующего вектора на каждом процессе
    double* y_total = new double[n];// вектор-результат
    double* a = new double[n * n];//исходная матрица
    if (rank == 0)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i * n + j] = rand() % 100;
            }
        }
        for (i = 0; i < n; i++)
        {
            x[i] = rand() % 100;
        }
    }
    double t = MPI_Wtime();
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a, n_partial * n, MPI_DOUBLE, a_partial, n_partial * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (i = 0; i < n_partial; i++)
    {
        for (j = 0; j < n; j++)
            y_partial[i] += a_partial[i * n + j] * x[j];
    }
    MPI_Gather(y_partial, n_partial, MPI_DOUBLE, y_total, n_partial, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    if (rank == 0)
    {
        std::cout << "time = " << t << "\n";
    }
    delete[] a_partial;
    delete[] a;
    delete[] x;
    delete[] y_partial;
    delete[] y_total;
}
void ParallelMatrix() //500х500 * 500х500
{
    int rank, size;
    int i, j, n = 500;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n_partial = n / size;
    double* a_partial = new double[n_partial * n];
    double* x = new double[n * n];
    double* y_partial = new double[n_partial * n];
    double* y_total = new double[n * n];
    double* a = new double[n * n];
    if (rank == 0)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i * n + j] = rand() % 100;
            }
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                x[i * n + j] = rand() % 100;
            }
        }
    }
    double t = MPI_Wtime();
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a, n_partial * n, MPI_DOUBLE, a_partial, n_partial * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (i = 0; i < n_partial; i++)
    {
        for (j = 0; j < n; j++)
            y_partial[i * n + j] += a_partial[i * n + j] * x[i * n + j];
    }
    MPI_Gather(y_partial, n_partial * n, MPI_DOUBLE, y_total, n_partial * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    if (rank == 0)
    {
        std::cout << "time = " << t << "\n";
    }
    delete[] a_partial;
    delete[] a;
    delete[] x;
    delete[] y_partial;
    delete[] y_total;
}
double rect_integral() { //Вычисление интеграла методом прямоугольников
    double b = 2.71828182845904;
    double a = 1;
    int n = 500;
    int ProcNumb;
    int ProcRank;
    double sum = 0;
    double summaGlobal;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNumb);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double t = MPI_Wtime();
    double s = (f(a) + f(b)) / 2;
    double h = (b - a) / n;
    for (int i = 0; i < n; ++i)
        s += f(a + i * h);
    sum = h * s;
    MPI_Reduce(&sum, &summaGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    std::cout << "time = " << t << " Значение = ";
    return sum;
}
double tr_integral() { //Вычисление интеграла методом трапеций
    double b = 2.71828182845904;
    double a = 1;
    int n = 500;
    double sum = 0;
    int ProcNumb;
    int ProcRank;
    double summaGlobal;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNumb);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double t = MPI_Wtime();
    const double width = (b - a) / n;
    for (int step = 0; step < n; step++) {
        const double x1 = a + step * width;
        const double x2 = a + (static_cast<__int64>(step) + 1) * width;
        sum += 0.5 * (x2 - x1) * (f(x1) + f(x2));
    }
    MPI_Reduce(&sum, &summaGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    std::cout << "time = " << t << " Значение = ";
    return sum;
}
void qsortRecursive(int* mas, int size) { //Сортировка
    int i = 0;
    int j = size - 1;
    int ProcNumb;
    int ProcRank;
    double summaGlobal;
    int* mass = new int[size]; int* c = new int[size]; int* cc = new int[size];
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNumb);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Scatter(mas, size, MPI_INT, mass, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int mid = mass[size / 2];
    do {
        while (mass[i] < mid) {
            i++;
        }
        while (mass[j] > mid) {
            j--;
        }
        if (i <= j) {
            int tmp = mass[i];
            mass[i] = mass[j];
            mass[j] = tmp;

            i++;
            j--;
        }
    } while (i <= j);
    if (j > 0) {
        qsortRecursive(mass, j + 1);
    }
    if (i < size) {
        qsortRecursive(&mass[i], size - i);
    }
    for (int i = 0; i < size; i++) {
        cc[i] = mass[i];
    };
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(cc, size, MPI_INT, c, size, MPI_INT, 0, MPI_COMM_WORLD);

}
int main(int argc, char** argv) {
    
    setlocale(LC_ALL, "RUS");
    printf("Быстрая сортировка 50000 элементов - ");
    double* pProcData;
    int     ProcDataSize;
    printf("Матрично-векторное умножение 500х500 * 500х1 - ");
    MPI_Init(&argc, &argv);
    ParallelMatrixV();
    printf("Матричное умножение 500х500 * 500х500 - ");
    ParallelMatrix();
    printf("Быстрая сортировка 50000 элементов - ");
    int S[5000];
    double t = MPI_Wtime();
    for (int j = 0; j < 5000; j++)
    {
        S[j] = rand() % 100;
    }
    qsortRecursive(S, 5000);
    t = MPI_Wtime() - t;
    std::cout << "time = " << t << "\n";
    printf("Метод прямоугольников - ");
    std::cout << rect_integral() << "\n";
    printf("Метод трапеций - ");
    std::cout << tr_integral() << "\n";
    getchar();
    MPI_Finalize();
}

//#include <iostream>
//#include <locale>
//#include <math.h>
//#include <stdio.h> 
//#include <stdlib.h> 
//#include <mpi.h>
//#include <vector>
//#include <random>
//
//using namespace std;
//
//#include "MPIUtilities.h"
//
////#define MPI_Matrix_multiply_vectors;
//
//void mpiMatrixMultiply(int size, int rank);
//vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int arows, int acols, int brows, int bcols, int procCount, int rank);
//vector<int> mpiSort(vector<int> a, int size, int rank);
//double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);
//double mpiRectangleMethod(double (*F)(double), double a, double b, int n);
//
//int MIN = 0;
//int MAX = 100;
//const double PI = 3.141592;
//
//// 0 -> Pi/6
//double F18(double x) {
//    return x * cos(3 * x);
//}
//
//#pragma region MPI_Matrix_multiply_variables
//const int arows = 500;
//const int acols = 500;
//const int brows = 500;
//const int bcols = 500;
//const int rowsPerProccess = 125; // rowsPerProccess = arows/size
//static int a[arows][acols];
//static int b[brows][bcols];
//int c[arows][bcols];
//#pragma endregion
//#pragma region MPI_Sort_variables
//const int sortNumsCount = 5000000;
//#pragma endregion
//#pragma region MPI_Math_variables
//
//#pragma endregion
//
//int main(int argc, char* argv[])
//{
//    double start, end; // Timer
//    int rank, size;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//#pragma region MPI_Matrix_multiply
//    if (rank == 0)
//    {
//        for (size_t i = 0; i < arows; i++)
//            for (size_t j = 0; j < acols; j++)
//                a[i][j] = randomInt();
//        for (size_t i = 0; i < brows; i++)
//            for (size_t j = 0; j < bcols; j++)
//                b[i][j] = randomInt();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(rank == 0)
//        start = MPI_Wtime();
//    mpiMatrixMultiply(size, rank);
//    if (rank == 0) {
//        end = MPI_Wtime();
//        cout << "Matrix mult "<< arows << "x" << acols << " * " << brows << "x" << bcols << " - " << abs(start - end) << " sec" << endl;
//    }
//#pragma endregion
//#pragma region MPI_Matrx_multiply_vectors
//#if defined MPI_Matrix_multiply_vectors
//    vector<vector<int>> a(8, vector<int>(8));
//    vector<vector<int>> b(8, vector<int>(8));
//    vector<vector<int>> c(8, vector<int>(8));
//
//    if (rank == 0)
//    {
//        a = vectorRandom(8, 8);
//        b = vectorRandom(8, 8);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (rank == 0)
//        start = MPI_Wtime();
//
//    c = mpiMatrixMultiply(a, b, 8, 8, 8, 8, size, rank);
//    if (rank == 0) {
//        end = MPI_Wtime();
//        cout << "Matrix mult " << 8 << "x" << 8 << " * " << 8 << "x" << 8 << " - " << abs(start - end) << " sec" << endl;
//    }
//#endif // MPI_Matrix_multiply_vectors
//#pragma endregion
//#pragma region MPI_Sort
//    static vector<int> t(sortNumsCount);
//    
//    if (rank == 0) {
//        t = vectorRandom(sortNumsCount);
//        start = MPI_Wtime();
//    }
//
//    auto c1 = mpiSort(t, size, rank);
//    
//    if (rank == 0) {
//        end = MPI_Wtime();
//        cout << "Sort " << sortNumsCount << " - " << abs(start - end) << " sec" << endl;
//    }
//#pragma endregion
//#pragma region MPI_Math
//
//#pragma endregion
//
//    MPI_Finalize();
//}
//vector<vector<int>> mpiMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b, int arows, int acols, int brows, int bcols, int procCount, int rank) {
//    if (arows % procCount != 0)
//        throw new exception("Количество строчек должно быть кратно количеству процессов");
//    if (acols != brows)
//        throw new exception("Нельзя перемножить");
//
//    int sum = 0;
//    int rowsPerProccess = arows / procCount;
//    int resultSize = arows * bcols / procCount;
//    
//    vector<vector<int>> aa(rowsPerProccess, vector<int>(acols));
//    vector<vector<int>> cc(rowsPerProccess, vector<int>(bcols));
//    vector<vector<int>> c(arows, vector<int>(bcols));
//    MPI_Scatter(&a.front(), rowsPerProccess * acols, MPI_INT, &aa.front(), rowsPerProccess * acols, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&b.front(), brows * bcols, MPI_INT, 0, MPI_COMM_WORLD);
//    
//    if (rank == 1)
//        cout << b[0][1] << endl;
//    
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i = 0; i < rowsPerProccess; i++)
//        for (int j = 0; j < bcols; j++) {
//            for (int s = 0; s < acols; s++)
//                sum += aa[i][s] * b[s][j];
//            cc[i][j] = sum;
//            sum = 0;
//        }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Gather(&cc.front(), resultSize, MPI_INT, &c.front(), resultSize, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//    return c;
//}
//void mpiMatrixMultiply(int size, int rank) {
//    if (arows % size != 0)
//        throw new exception("Количество строчек должно быть кратно количеству процессов");
//    if (acols != brows)
//        throw new exception("Нельзя перемножить");
//
//    int sum = 0;
//    int matrixSize = acols;    
//    int partSize = arows * rowsPerProccess;
//    int resultSize = arows * bcols / size;
//    auto aa = new int[rowsPerProccess][acols];
//    auto cc = new int[rowsPerProccess][bcols];
//    
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Scatter(a, partSize, MPI_INT, aa, partSize, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Bcast(b, brows * bcols, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    for (int i = 0; i < rowsPerProccess; i++)
//        for (int j = 0; j < bcols; j++) {
//            for (int s = 0; s < matrixSize; s++)
//                sum += aa[i][s] * b[s][j];
//            cc[i][j] = sum;
//            sum = 0;
//        }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    MPI_Gather(cc, resultSize, MPI_INT, c, resultSize, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//}
//vector<int> mpiSort(vector<int> a, int size, int rank)
//{
//    auto numsPerProccess = a.size() / size;
//    static vector<int> c(a.size());
//    static vector<int> aa(numsPerProccess);
//    MPI_Barrier(MPI_COMM_WORLD);
//    
//    MPI_Scatter(&a.front(), numsPerProccess, MPI_INT, &aa.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//    
//    sort(aa.begin(), aa.end());
//    MPI_Barrier(MPI_COMM_WORLD);
//    
//    MPI_Gather(&aa.front(), numsPerProccess, MPI_INT, &c.front(), numsPerProccess, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rank == 0)
//    {
//        sort(c.begin(), c.end());
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    return c;
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