#include <iostream>
#include <locale>
#include <math.h>
#include <random>
#include <iomanip>
#include <Windows.h>
#include <stdio.h> 
#include <time.h> 
#include <vector>
#include <format>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

#include "Openmp.h"
#include "Serial.h"
#include "PPL.h"
#include "Utilities.h"

void matrixMultComparsResult();
void mathCompars(int stepCount=1000);
void sortCompars(int vectorLenght=30);
void matrixMultCompars(int matrixSize=500);
void PPLTest();
void OpenMpTest();
void serialTest();
void getAccuracy(int n);

const double PI = 3.141592;
int MIN = 0;
int MAX = 100;

double F(double x) {
    return(1/log(x));
}
// 0 -> Pi/6
double F18(double x) {
    return x * cos(3 * x);
}
// 0 -> Pi
double F19(double x) {
    return x * x * sin(x);
}
template<typename Func, typename... Args>
void printExecutionTime(string msg, Func f, Args&&... args) {
    clock_t start = clock();
    f(args...);
    clock_t end = clock();
    double seconds = (double)((end - start) * 1000) / CLOCKS_PER_SEC;
    std::cout << "Execution time: " << seconds << "ms " << msg << std::endl;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << setprecision(16);

    getAccuracy(10);

    //serialTest();
    //PPLTest();
    //OpenMpTest();
}
// n - количество знаков после запятой
void getAccuracy(int n) {
    double prevR = -10, R = 10;
    for (int i = 10; abs(R - prevR) > pow(10, -n); i += 1)
    {
        prevR = R;
        R = rectangleMethod(F18, 0, PI / 6, i);
        cout << "rectangleMethod при " << i << " итерациях, расхождение - "<< abs(R - prevR) << endl;
    }
}
void mathCompars(int stepCount) {
    cout << "rectangleMethod " << rectangleMethod(F18, 0, PI / 6, stepCount) << endl;
    cout << "trapeziodalMethod " << trapeziodalMethod(F18, 0, PI / 6, stepCount) << endl;
    
    printExecutionTime("Метод прямоугольников", rectangleMethod, F18, 0, PI / 6, stepCount);
    printExecutionTime("Метод трапеций", trapeziodalMethod, F18, 0, PI / 6, stepCount);
    cout << endl;

    cout << "parallelRectangleMethod " << parallelRectangleMethod(F18, 0, PI / 6, stepCount) << endl;
    cout << "parallelTrapeziodalMethod " << parallelTrapeziodalMethod(F18, 0, PI / 6, stepCount) << endl;

    printExecutionTime("Параллельный метод прямоугольников", parallelRectangleMethod, F18, 0, PI / 6, stepCount);
    printExecutionTime("Параллельный метод трапеций", parallelTrapeziodalMethod, F18, 0, PI / 6, stepCount);
    cout << endl;

    cout << "openmpRectangleMethod " << openmpRectangleMethod(F18, 0, PI / 6, stepCount) << endl;
    cout << "openmpTrapeziodalMethod " << openmpTrapeziodalMethod(F18, 0, PI / 6, stepCount) << endl;

    printExecutionTime("OpenMP метод прямоугольников", openmpRectangleMethod, F18, 0, PI / 6, stepCount);
    printExecutionTime("OpenMP метод трапеций", openmpTrapeziodalMethod, F18, 0, PI / 6, stepCount);
    cout << endl;
}
void sortCompars(int vectorLenght) {
    auto vec = vectorRandom(vectorLenght);
    auto temp = vec;

    printVector(temp);
    printExecutionTime("Сортировка вставками", sortInsertion, temp);
    printVector(temp);

    temp = vec;
    printVector(temp);
    printExecutionTime("Системная сортировка", sortSystem, temp);
    printVector(temp);

    temp = vec;
    printVector(temp);
    printExecutionTime("Параллельная сортировка вставками", parallelSortInsertion, temp);
    printVector(temp);

    temp = vec;
    printVector(temp);
    printExecutionTime("Параллельная сортировка вставками", openmpSortInsertion, temp);
    printVector(temp);
}
void matrixMultCompars(int matrixSize) {
    auto matrixSizeString = std::to_string(matrixSize);
    
    auto matrixVectorMult = matrixSizeString + "x" + matrixSizeString + " * " + matrixSizeString + "x1";
    auto matrixMatrixMult = matrixSizeString + "x" + matrixSizeString + " * " + matrixSizeString + "x" + matrixSizeString;
    
    printExecutionTime(matrixVectorMult, matrixSimpleMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, 1));
    printExecutionTime(matrixMatrixMult, matrixSimpleMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, matrixSize));

    printExecutionTime("PPL " + matrixVectorMult, parallelMatrixMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, 1));
    printExecutionTime("PPL " + matrixMatrixMult, parallelMatrixMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, matrixSize));

    printExecutionTime("OpenMP " + matrixVectorMult, openmpMatrixMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, 1));
    printExecutionTime("OpenMP " + matrixMatrixMult, openmpMatrixMultiply, vectorRandom(matrixSize, matrixSize), vectorRandom(matrixSize, matrixSize));
}
void matrixMultComparsResult() {
    vector<vector<int>> a {
        { 1, 2, 3, 4, 5 },
        { 2, 2, 2, 2, 2 },
        { 3, 1, 2, 4, 2 },
        { -2, -2, -2, -2, -2 },
        { 5, 4, 3, 2, 1 }
    };

    vector<vector<int>> b {
        { 1, 3, 3, 4, 5 },
        { 1, 3, 1, 1, 1 },
        { 3, 3, 2, 4, 2 },
        { -2, -3, -2, -2, -1 },
        { 5, 3, 3, 2, 1 }
    };

    auto c = matrixSimpleMultiply(a, b);
    cout << "Default" << endl;
    printVector(c);

    c = parallelMatrixMultiply(a, b);
    cout << "PPL" << endl;
    printVector(c);

    c = openmpMatrixMultiply(a, b);
    cout << "OpenMp" << endl;
    printVector(c);
}
void OpenMpTest() {
    cout << "OpenMP\r\n";

    printExecutionTime("500x500 * 500x1", openmpMatrixMultiply, vectorRandom(500, 500), vectorRandom(500, 1));
    printExecutionTime("500x500 * 500x500", openmpMatrixMultiply, vectorRandom(500, 500), vectorRandom(500, 500));

    auto temp = vectorRandom(10000);
    printExecutionTime("Параллельная сортировка вставками", openmpSortInsertion, temp);

    cout << "parallelRectangleMethod " << openmpRectangleMethod(F18, 0, PI / 6, 100) << endl;
    cout << "parallelTrapeziodalMethod " << openmpTrapeziodalMethod(F18, 0, PI / 6, 100) << endl;

    printExecutionTime("Параллельный метод прямоугольников", openmpRectangleMethod, F18, 0, PI / 6, 1000);
    printExecutionTime("Параллельный метод трапеций", openmpTrapeziodalMethod, F18, 0, PI / 6, 1000);

    cout << endl;
}
void PPLTest() {
    cout << "PPL\r\n";

    printExecutionTime("500x500 * 500x1", parallelMatrixMultiply, vectorRandom(500, 500), vectorRandom(500, 1));
    printExecutionTime("500x500 * 500x500", parallelMatrixMultiply, vectorRandom(500, 500), vectorRandom(500, 500));

    auto temp = vectorRandom(10000);
    printExecutionTime("Параллельная сортировка вставками", parallelSortInsertion, temp);

    cout << "parallelRectangleMethod " << parallelRectangleMethod(F18, 0, PI / 6, 100) << endl;
    cout << "parallelTrapeziodalMethod " << parallelTrapeziodalMethod(F18, 0, PI / 6, 100) << endl;
    
    printExecutionTime("Параллельный метод прямоугольников", parallelRectangleMethod, F18, 0, PI / 6, 1000);
    printExecutionTime("Параллельный метод трапеций", parallelTrapeziodalMethod, F18, 0, PI / 6, 1000);

    cout << endl;
}
void serialTest() {
    cout << "Serial\r\n";

    printExecutionTime("500x500 * 500x1", matrixSimpleMultiply, vectorRandom(500, 500), vectorRandom(500, 1));
    printExecutionTime("500x500 * 500x500", matrixSimpleMultiply, vectorRandom(500, 500), vectorRandom(500, 500));

    auto vec = vectorRandom(10000);
    auto temp = vec;
    printExecutionTime("Сортировка вставками", sortInsertion, temp);

    temp = vec;
    printExecutionTime("Системная сортировка", sortSystem, temp);

    cout << "trapeziodalMethod " << trapeziodalMethod(F18, 0, PI / 6, 100) << endl;
    cout << "rectangleMethod " << rectangleMethod(F18, 0, PI / 6, 100) << endl;

    printExecutionTime("Метод прямоугольников", rectangleMethod, F18, 0, PI / 6, 1000);
    printExecutionTime("Метод трапеций", trapeziodalMethod, F18, 0, PI / 6, 1000);

    cout << endl;
}