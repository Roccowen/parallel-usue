#include <vector>
#include <random>
#include <ppl.h>
#include "PPL.h"

using namespace std;

void parallelSortInsertion(vector<int>& arr)
{
    size_t length = arr.size();
    int temp;
    bool sorted = false;

    while (!sorted)
    {
        sorted = true;
        Concurrency::parallel_for(size_t(0), length - 1, [&](size_t i)
            {
                if (arr[i] > arr[i + 1])
                {
                    temp = arr[i + 1];
                    arr[i + 1] = arr[i];
                    arr[i] = temp;
                    sorted = false;
                }
            });
    }
}
std::vector<std::vector<int>> parallelMatrixMultiply(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b)
{
    if (a[0].size() != b.size())
        throw new exception("error");

    size_t rows = a.size();
    size_t cols = b[0].size();
    size_t size = a[0].size();
    vector<vector<int>> c(rows, vector<int>(cols));

    Concurrency::parallel_for(size_t(0), rows, [&](size_t i)
        {
            for (size_t j = 0; j < cols; j++)
            {
                double temp = 0;
                for (int k = 0; k < size; k++)
                {
                    temp += a[i][k] * b[k][j];
                }
                c[i][j] = temp;
            }
        });
    return c;
}
double parallelRectangleMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double start = a + step;
    double stop = b;
    std::atomic<double> result = 0;

    Concurrency::parallel_for(0, n, 1, [&](int step_num)
        {
            result += F(step_num * step + start);
        }
    );
        
    return result * step;
}
double parallelTrapeziodalMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double start = a + step;
    double stop = b;
    std::atomic<double> result = 0;

    Concurrency::parallel_for(0, n, 1, [&](int step_num)
        {
            result += F(step_num * step + start);
        }
    );

    return (result + F(a) / 2 + F(b) / 2) * step;
}