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
extern "C" {
#include <cblas.h>
}
#endif

namespace amici {

void amici_dgemm(
    BLASLayout layout, BLASTranspose TransA, BLASTranspose TransB, int const M,
    int const N, int const K, double const alpha, double const* A,
    int const lda, double const* B, int const ldb, double const beta, double* C,
    int const ldc
) {
    cblas_dgemm(
        (CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
        M, N, K, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

void amici_dgemv(
    BLASLayout layout, BLASTranspose TransA, int const M, int const N,
    double const alpha, double const* A, int const lda, double const* X,
    int const incX, double const beta, double* Y, int const incY
) {
    cblas_dgemv(
        (CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A, lda, X,
        incX, beta, Y, incY
    );
}

void amici_daxpy(
    int n, double alpha, double const* x, int const incx, double* y, int incy
) {
    cblas_daxpy(n, alpha, x, incx, y, incy);
}

} // namespace amici
