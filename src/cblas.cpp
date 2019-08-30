/**
 * @file   cblas.cpp
 * @brief  BLAS routines required by AMICI
 *
 **/

#include "amici/cblas.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif defined(AMICI_BLAS_MKL)
#include <mkl.h>
#else
extern "C"
{
   #include <cblas.h>
}
#endif

namespace amici {

void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    cblas_dgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA,
                (CBLAS_TRANSPOSE)TransB, M, N, K, alpha, A, lda, B, ldb, beta,
                C, ldc);
}

void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    cblas_dgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A,
                lda, X, incX, beta, Y, incY);
}

void amici_daxpy(int n, double alpha, const double *x, const int incx, double *y, int incy) {
    cblas_daxpy(n, alpha, x, incx, y, incy);
}

} // namespace amici
