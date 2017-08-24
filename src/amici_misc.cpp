#include "include/amici_misc.h"

#include <cstdio>
#include <cstring>

void fillArray(double *destination, int count, double value) {
    for (int i = 0; i < count; ++i)
        destination[i] = value;
}

double sum(double const *array, int numElements) {
    double sum = 0;
    for (int i = 0; i < numElements; ++i) {
        sum += array[i];
    }
    return sum;
}

void zeros(double *destination, int count) {
    memset(destination, 0, sizeof(double) * count);
}

void ones(double *destination, int count) { fillArray(destination, count, 1); }

void linSpace(double *destination, double from, double to, int numValues) {
    double delta = (to - from) / (numValues - 1);
    int i;
    for (i = 0; i < numValues; ++i) {
        destination[i] = from + i * delta;
    }
}

double *linSpaceAlloc(double from, double to, int numValues) {
    double *destination = new double[numValues];
    linSpace(destination, from, to, numValues);
    return destination;
}

void printArray(double const *array, int numElements) {
    printfArray(array, numElements, "%e\t");
}

void printfArray(double const *array, int numElements, char const *format) {
    for (int i = 0; i < numElements; ++i) {
        printf(format, array[i]);
    }
}
