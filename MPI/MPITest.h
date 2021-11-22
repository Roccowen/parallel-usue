#ifndef MPITEST

#define MPITEST
#include <vector>

/// <summary>
/// Произведение матриц
/// </summary>
/// <param name="a">Первая матрица</param>
/// <param name="b">Вторая матрица</param>
/// <param name="aRows">Количество строчек в первой матрице</param>
/// <param name="aCols">Количество столбцов в первой матрице</param>
/// <param name="bRows">Количество строчек во второй матрице</param>
/// <param name="bCols">Количество столбцов во второй матрице</param>
/// <param name="size">Количество потоков</param>
/// <returns></returns>
int** mpiMatrixMultiply(int** a, int** b, int aRows, int aCols, int bRows, int bCols, int size);
void mpiMatrixMultiplyTest(int size);
void mpiSort(std::vector<int>& arr);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !SERIAL
