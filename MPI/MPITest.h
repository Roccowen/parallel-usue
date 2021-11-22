#ifndef MPITEST

#define MPITEST
#include <vector>

/// <summary>
/// ������������ ������
/// </summary>
/// <param name="a">������ �������</param>
/// <param name="b">������ �������</param>
/// <param name="aRows">���������� ������� � ������ �������</param>
/// <param name="aCols">���������� �������� � ������ �������</param>
/// <param name="bRows">���������� ������� �� ������ �������</param>
/// <param name="bCols">���������� �������� �� ������ �������</param>
/// <param name="size">���������� �������</param>
/// <returns></returns>
int** mpiMatrixMultiply(int** a, int** b, int aRows, int aCols, int bRows, int bCols, int size);
void mpiMatrixMultiplyTest(int size);
void mpiSort(std::vector<int>& arr);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !SERIAL
