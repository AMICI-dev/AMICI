#include "amici/misc.h"
#include "amici/amici.h"
#include "amici/symbolic_functions.h"

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

/** Checks the values in an array for NaNs and Infs
 *
 * @param N number of elements in array
 * @param array array
 * @param fun name of calling function
 * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found, AMICI_SUCCESS otherwise
 */
int checkFinite(const int N,const realtype *array, const char* fun){
    for(int idx = 0; idx < N; idx++) {
        if(isNaN(array[idx])) {
            warnMsgIdAndTxt("AMICI:mex:NaN","AMICI encountered a NaN value at index %i of %i in %s! Trying to recover ... ",idx,N,fun);
            return(AMICI_RECOVERABLE_ERROR);
        }
        if(isInf(array[idx])) {
            warnMsgIdAndTxt("AMICI:mex:Inf","AMICI encountered an Inf value at index %i of %i in %s! Trying to recover ... ",idx,N,fun);
            return(AMICI_RECOVERABLE_ERROR);
        }
    }
    return(AMICI_SUCCESS);
}

} // namespace amici
