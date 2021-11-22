#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <Windows.h>
#include "MPIUtilities.h"

using namespace std;

extern int MIN;
extern int MAX;

int randomInt() {
    static std::random_device rd;                            // only used once to initialise (seed) engine
    static std::mt19937 rng(rd());                           // random-number engine used (Mersenne-Twister in this case)
    static std::uniform_int_distribution<int> uni(MIN, MAX); // guaranteed unbiased

    return uni(rng);
}
void printArray(int* arr, int arrlenght) {
    for (size_t i = 0; i < arrlenght; i++)
        cout << arr[i] << " ";
    cout << endl;
}
void printArray(int** arr, int arrRows, int arrColumns) {
    for (size_t i = 0; i < arrRows; i++ ) {
        for (size_t j = 0; j < arrColumns; j++ )
            cout << arr[i][j] << " ";
        cout << endl;
    }
}
int* arrayRandom(int count)
{
    auto arr = new int[count];
    for (size_t i = 0; i < count; i++)
        arr[i] = randomInt();

    return arr;
}
int** arrayRandom(int rows, int cols)
{
    int** arr2d = new int* [rows];
    for (int i = 0; i < rows; ++i)
        arr2d[i] = new int[cols];

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            arr2d[i][j] = randomInt();

    return arr2d;
}