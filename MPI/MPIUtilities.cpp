#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <Windows.h>
#include "MPIUtilities.h"

using namespace std;

//extern int MIN;
//extern int MAX;

int randomInt() {
    static std::random_device rd;                            // only used once to initialise (seed) engine
    static std::mt19937 rng(rd());                           // random-number engine used (Mersenne-Twister in this case)
    //static std::uniform_int_distribution<int> uni(MIN, MAX); // guaranteed unbiased
    static std::uniform_int_distribution<int> uni(1, 100); // guaranteed unbiased

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
void printVector(vector<int> arr) {
    for (int i : arr)
        cout << i << " ";
    cout << endl;
}
void printVector(std::vector<std::vector<int>> matrix) {
    for (auto i : matrix) {
        for (auto j : i)
            cout << j << " ";
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
vector<int> vectorRandom(int count)
{
    vector<int> random_vector(count);
    for (size_t i = 0; i < count; i++)
        random_vector[i] = randomInt();

    return random_vector;
}
vector<vector<int>> vectorRandom(int rows, int cols)
{
    vector<vector<int> > random_vector(rows, vector<int>(cols));
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            random_vector[i][j] = randomInt();

    return random_vector;
}