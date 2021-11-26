#ifndef MPIUTILITIES

#define MPIUTILITIES
#include <vector>


int randomInt();
void printArray(int* arr, int arrlenght);
void printArray(int** arr, int arrRows, int arrColumns);
void printVector(std::vector<int> arr);
int* arrayRandom(int count);
int** arrayRandom(int rows, int cols);
std::vector<int> vectorRandom(int count);
std::vector<std::vector<int>> vectorRandom(int rows, int cols);

#endif // !MPIUTILITIES