#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void ones(double *destination, int count);
EXTERNC void zeros(double *destination, int count);
EXTERNC void fillArray(double *destination, int count, double value);
EXTERNC double sum(double const *array, int numElements);
EXTERNC void linSpace(double *destination, double from, double to, int numValues);
EXTERNC double *linSpaceAlloc(double from, double to, int numValues);
EXTERNC void printArray(double const *array, int numElements);
EXTERNC void printfArray(double const *array, int numElements, char const *format);

#endif // AMICI_MISC_H
