#ifndef MPITEST

#define MPITEST
#include <vector>

std::vector<std::vector<int>> mpiMatrixMultiply(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b);
void mpiSort(std::vector<int>& arr);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !SERIAL
