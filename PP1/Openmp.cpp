#include <vector>
#include <random>
#include <omp.h> 
#include <thread>
#include <vector>

#include "Openmp.h"

using namespace std;

vector<vector<int>> openmpMatrixMultiply(vector<vector<int>> a, vector<vector<int>> b) {
	if (a[0].size() != b.size())
		throw new exception("error");

	int rows = a.size();
	int cols = b[0].size();
	int size = a[0].size();
	vector<vector<int>> c(rows, vector<int>(cols));
	int i, j, s;
    int temp;

#pragma omp parallel private(i, j, s)
    {
#pragma omp for 
        for (i = 0; i < rows; i++)
            for (j = 0; j < cols; j++) {
                c[i][j] = 0;
                for (s = 0; s < size; s++)
                    c[i][j] += a[i][s] * b[s][j];
            }
	}
	return c;
}

void openmpSortInsertion(vector<int>& arr)
{
    size_t length = arr.size();
    int temp = 0;
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;
#pragma omp parallel shared(arr, length)
#pragma omp for    
        for (int i = 0; i < length - 1; i++)
        {
            if (arr[i] > arr[i + 1])
            {
#pragma omp critical
                {
                temp = arr[i + 1];
                arr[i + 1] = arr[i];
                arr[i] = temp;
                sorted = false;
                }
            }
        }
    }
}

double openmpRectangleMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double start = a + step;
    double stop = b;
    double result = 0;

#pragma omp parallel reduction(+: result) shared(step, start, n)
#pragma omp for
    for (int step_num = 0; step_num < n; step_num++)
        result += F(step_num * step + start);

    return result * step;
}
double openmpTrapeziodalMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double start = a + step;
    double stop = b;
    double result = 0;

#pragma omp parallel reduction(+: result) shared(step, start, n)
#pragma omp for
    for (int step_num = 0; step_num < n; step_num++)
        result += F(step_num * step + start);
    
    return (result + F(a) / 2 + F(b) / 2) * step;
}