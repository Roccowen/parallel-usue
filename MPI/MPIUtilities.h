#ifndef MPIUTILITIES

#define MPIUTILITIES

int randomInt();
void printArray(int* arr, int arrlenght);
void printArray(int** arr, int arrRows, int arrColumns);
int* arrayRandom(int count);
int** arrayRandom(int rows, int cols);

#endif // !MPIUTILITIES