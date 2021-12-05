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
void PPLTest(int a1rows = 500, int a1cols = 500, int b1rows = 500, int b1cols = 1, int a2rows = 500, int a2cols = 500, int b2rows = 500, int b2cols = 500, int sortCount = 10000, int mathStep = 1000, int repeatCount = 1);
void OpenMpTest(int a1rows = 500, int a1cols = 500, int b1rows = 500, int b1cols = 1, int a2rows = 500, int a2cols = 500, int b2rows = 500, int b2cols = 500, int sortCount = 10000, int mathStep = 1000, int repeatCount = 1);
void serialTest(int a1rows = 500, int a1cols = 500, int b1rows = 500, int b1cols = 1, int a2rows = 500, int a2cols = 500, int b2rows = 500, int b2cols = 500, int sortCount = 10000, int mathStep = 1000, int repeatCount = 1);
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
    //PPLTest(1800, 1800, 1800, 1, 1800, 1800, 1800, 1800, 1800, 10000, 3);
    //PPLTest(2800, 2800, 2800, 1, 2800, 2800, 2800, 2800, 11800, 20000, 3);
    //PPLTest(3800, 3800, 3800, 1, 3800, 3800, 3800, 3800, 21800, 30000, 3);
    //PPLTest(4800, 4800, 4800, 1, 4800, 4800, 4800, 4800, 31800, 40000, 3);
    //PPLTest(5800, 5800, 5800, 1, 5800, 5800, 5800, 5800, 41800, 50000, 3);
    //OpenMpTest(1800, 1800, 1800, 1, 1800, 1800, 1800, 1800, 1800, 10000, 3);
    //OpenMpTest(2800, 2800, 2800, 1, 2800, 2800, 2800, 2800, 11800, 20000, 3);
    //OpenMpTest(3800, 3800, 3800, 1, 3800, 3800, 3800, 3800, 21800, 30000, 3);
    //OpenMpTest(4800, 4800, 4800, 1, 4800, 4800, 4800, 4800, 31800, 40000, 3);
    //OpenMpTest(5800, 5800, 5800, 1, 5800, 5800, 5800, 5800, 41800, 50000, 3);
    serialTest(1800, 1800, 1800, 1, 1800, 1800, 1800, 1800, 1800, 10000, 1);
    serialTest(2800, 2800, 2800, 1, 2800, 2800, 2800, 2800, 11800, 20000, 1);
    serialTest(3800, 3800, 3800, 1, 3800, 3800, 3800, 3800, 21800, 30000, 1);
    serialTest(4800, 4800, 4800, 1, 4800, 4800, 4800, 4800, 31800, 40000, 1);
    serialTest(5800, 5800, 5800, 1, 5800, 5800, 5800, 5800, 41800, 50000, 1);

    //getAccuracy(10);

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

void OpenMpTest(int a1rows, int a1cols, int b1rows, int b1cols, int a2rows, int a2cols, int b2rows, int b2cols, int sortCount, int mathStep, int repeatCount) {
    cout << "OpenMP\r\n";

    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult1 = std::to_string(a1rows) + "x" + std::to_string(a1cols) + " * " + std::to_string(b1rows) + "x" + std::to_string(b1cols);
        printExecutionTime(matrixMult1, openmpMatrixMultiply, vectorRandom(a1rows, a1cols), vectorRandom(b1rows, b1cols));
    }
    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult2 = std::to_string(a2rows) + "x" + std::to_string(a2cols) + " * " + std::to_string(b2rows) + "x" + std::to_string(b2cols);
        printExecutionTime(matrixMult2, openmpMatrixMultiply, vectorRandom(a2rows, a2cols), vectorRandom(b2rows, b2cols));
    }
    auto vec = vectorRandom(sortCount);
    auto temp = vec;
    for (size_t i = 0; i < repeatCount; i++) {
        temp = vec;
        printExecutionTime("Параллельная сортировка вставками", openmpSortInsertion, temp);
    }
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Параллельный метод прямоугольников", openmpRectangleMethod, F18, 0, PI / 6, mathStep);
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Параллельный метод трапеций", openmpTrapeziodalMethod, F18, 0, PI / 6, mathStep);
    cout << endl;
}
void PPLTest(int a1rows, int a1cols, int b1rows, int b1cols, int a2rows, int a2cols, int b2rows, int b2cols, int sortCount, int mathStep, int repeatCount) {
    cout << "PPL\r\n";

    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult1 = std::to_string(a1rows) + "x" + std::to_string(a1cols) + " * " + std::to_string(b1rows) + "x" + std::to_string(b1cols);
        printExecutionTime(matrixMult1, parallelMatrixMultiply, vectorRandom(a1rows, a1cols), vectorRandom(b1rows, b1cols));
    }
    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult2 = std::to_string(a2rows) + "x" + std::to_string(a2cols) + " * " + std::to_string(b2rows) + "x" + std::to_string(b2cols);
        printExecutionTime(matrixMult2, parallelMatrixMultiply, vectorRandom(a2rows, a2cols), vectorRandom(b2rows, b2cols));
    }
    auto vec = vectorRandom(sortCount);
    auto temp = vec;
    for (size_t i = 0; i < repeatCount; i++) {
        temp = vec;
        printExecutionTime("Параллельная сортировка вставками", parallelSortInsertion, temp);
    }
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Параллельный метод прямоугольников", parallelRectangleMethod, F18, 0, PI / 6, mathStep);
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Параллельный метод трапеций", parallelTrapeziodalMethod, F18, 0, PI / 6, mathStep);
    cout << endl;
}
void serialTest(int a1rows, int a1cols, int b1rows, int b1cols, int a2rows, int a2cols, int b2rows, int b2cols, int sortCount, int mathStep, int repeatCount) {
    cout << "Serial\r\n";

    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult1 = std::to_string(a1rows) + "x" + std::to_string(a1cols) + " * " + std::to_string(b1rows) + "x" + std::to_string(b1cols);
        printExecutionTime(matrixMult1, matrixSimpleMultiply, vectorRandom(a1rows, a1cols), vectorRandom(b1rows, b1cols));
    }
    for (size_t i = 0; i < repeatCount; i++)
    {
        auto matrixMult2 = std::to_string(a2rows) + "x" + std::to_string(a2cols) + " * " + std::to_string(b2rows) + "x" + std::to_string(b2cols);
        printExecutionTime(matrixMult2, matrixSimpleMultiply, vectorRandom(a2rows, a2cols), vectorRandom(b2rows, b2cols));
    }
    auto vec = vectorRandom(sortCount);
    auto temp = vec;
    for (size_t i = 0; i < repeatCount; i++) {
        temp = vec;
        printExecutionTime("Сортировка вставками", sortInsertion, temp);
    }    
    for (size_t i = 0; i < repeatCount; i++)
    {
        temp = vec;
        printExecutionTime("Системная сортировка", sortSystem, temp);
    }
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Метод прямоугольников", rectangleMethod, F18, 0, PI / 6, mathStep);
    for (size_t i = 0; i < repeatCount; i++)
        printExecutionTime("Метод трапеций", trapeziodalMethod, F18, 0, PI / 6, mathStep);
    cout << endl;
}