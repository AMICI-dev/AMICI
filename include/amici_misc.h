#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include<algorithm>
void zeros(double *destination, int count);
void fillArray(double *destination, int count, double value);
void printArray(double const *array, int numElements);
void printfArray(double const *array, int numElements,
                         char const *format);

#endif // AMICI_MISC_H
