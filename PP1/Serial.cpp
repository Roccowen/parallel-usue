#include <vector>
#include <random>
#include "Serial.h"

using namespace std;

vector<vector<int>> matrixSimpleMultiply(vector<vector<int>> a, vector<vector<int>> b)
{
    if (a[0].size() != b.size())
        throw new exception("error");

    int rows = a.size();
    int cols = b[0].size();
    int size = a[0].size();
    vector<vector<int> > c(rows, vector<int>(cols));

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            for (size_t s = 0; s < size; s++)
                c[i][j] += a[i][s] * b[s][j];

    return c;
}
void sortReplace(vector<int>& arr) {
    int buff;
    int arr_size = arr.size();
    for (size_t i = 1; i < arr_size; i++)
    {
        buff = arr[i];
        for (size_t j = i - 1; j >= 0 && arr[j] > buff; j--)
        {
            arr[j + 1] = arr[j];
            arr[j + 1] = buff;
        }
    }
}
void sortSystem(vector<int>& arr)
{
    sort(arr.begin(), arr.end());
}
void sortInsertion(vector<int>& arr)
{
    int length = arr.size();
    int* arrayPtr = &arr[0];
    int temp,
        item;
    for (int counter = 1; counter < length; counter++)
    {
        temp = arrayPtr[counter];
        item = counter - 1;
        while (item >= 0 && arrayPtr[item] > temp)
        {
            arrayPtr[item + 1] = arrayPtr[item];
            arrayPtr[item] = temp;
            item--;
        }
    }
}
double trapeziodalMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double result = 0;
    for (double x = a + step; x < b; x += step)
        result += F(x);
    return (result + F(a) / 2 + F(b) / 2) * step;
}
double rectangleMethod(double (*F)(double), double a, double b, int n) {
    double step = (b - a) / n;
    double result = 0;
    for (double x = a + step; x <= b; x += step)
        result += F(x);
    return result * step;
}