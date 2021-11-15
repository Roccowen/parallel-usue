#ifndef PPL

#define PPL
#include <vector>

void parallelSortInsertion(std::vector<int>& arr);
std::vector<std::vector<int>> parallelMatrixMultiply(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b);
double parallelRectangleMethod(double (*F)(double), double a, double b, int n);
double parallelTrapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !PPL