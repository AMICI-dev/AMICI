#ifndef AMICI_CBLAS_H
#define AMICI_CBLAS_H

#include "amici/defines.h"

namespace amici {

void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 int M, int N, double alpha, const double *A,
                 int lda, const double *X, int incX,
                 double beta, double *Y, int incY);

void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, int M, int N,
                 int K, double alpha, const double *A,
                 int lda, const double *B, int ldb,
                 double beta, double *C, int ldc);

void amici_daxpy(int n, double alpha, const double *x, int incx, double *y, int incy);

} // namespace amici

#endif // AMICI_CBLAS_H
