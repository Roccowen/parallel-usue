#ifndef SERIAL

#define SERIAL
#include <vector>

std::vector<std::vector<int>> matrixSimpleMultiply(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b);
void sortReplace(std::vector<int>& arr);
void sortSystem(std::vector<int>& arr);
void sortInsertion(std::vector<int>& arr);
double rectangleMethod(double (*F)(double), double a, double b, int n);
double trapeziodalMethod(double (*F)(double), double a, double b, int n);

#endif // !SERIAL