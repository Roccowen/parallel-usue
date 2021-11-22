#ifndef MPITEST

#define MPITEST
#include <vector>



void mpiSort(std::vector<int>& arr);
double mpiRectangleMethod(double (*F)(double), double a, double b, int n);
double mpiTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !SERIAL
