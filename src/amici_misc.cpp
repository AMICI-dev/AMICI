#include "include/amici_misc.h"

#include <cstdio>
#include <cstring>

namespace amici {

void fillArray(double *destination, int count, double value) {
    for (int i = 0; i < count; ++i)
        destination[i] = value;
}

void zeros(double *destination, int count) {
    memset(destination, 0, sizeof(double) * count);
}

void printArray(double const *array, int numElements) {
    printfArray(array, numElements, "%e\t");
}

void printfArray(double const *array, int numElements, char const *format) {
    for (int i = 0; i < numElements; ++i) {
        printf(format, array[i]);
    }
}

} // namespace amici
