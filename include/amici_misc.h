#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void zeros(double *destination, int count);
EXTERNC void fillArray(double *destination, int count, double value);
EXTERNC void printArray(double const *array, int numElements);
EXTERNC void printfArray(double const *array, int numElements,
                         char const *format);

#endif // AMICI_MISC_H
