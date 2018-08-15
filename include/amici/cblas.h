#ifndef AMICI_CBLAS_H
#define AMICI_CBLAS_H

#include "amici/defines.h"

namespace amici {

void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);

void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

void amici_daxpy(int n, double alpha, const double *x, const int incx, double *y, int incy);

} // namespace amici

#endif // AMICI_CBLAS_H
