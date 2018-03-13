#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include "amici/defines.h"

#include<algorithm>


namespace amici {

void zeros(double *destination, int count);
void fillArray(double *destination, int count, double value);
void printArray(double const *array, int numElements);
void printfArray(double const *array, int numElements,
                         char const *format);
int isFinite(const int N,const realtype *array, const char* fun);

} // namespace amici
#endif // AMICI_MISC_H
