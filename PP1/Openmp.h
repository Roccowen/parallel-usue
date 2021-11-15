#ifndef Openmp

#define Openmp
#include <vector>

std::vector<std::vector<int>> openmpMatrixMultiply(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b);
void openmpSortInsertion(std::vector<int>& arr);
double openmpRectangleMethod(double (*F)(double), double a, double b, int n);
double openmpTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !PPL