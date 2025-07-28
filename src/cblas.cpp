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

#ifndef BLAS_PREFIX
#define BLAS_PREFIX cblas_
#endif

#ifndef BLAS_SUFFIX
#define BLAS_SUFFIX
#endif

#define BLAS_CONCAT2(a, b, c) a##b##c
#define BLAS_CONCAT(a, b, c) BLAS_CONCAT2(a, b, c)
#define BLAS_FUNC(name) BLAS_CONCAT(BLAS_PREFIX, name, BLAS_SUFFIX)

namespace amici {

void amici_dgemm(
    BLASLayout layout, BLASTranspose TransA, BLASTranspose TransB, int const M,
    int const N, int const K, double const alpha, double const* A,
    int const lda, double const* B, int const ldb, double const beta, double* C,
    int const ldc
) {

    BLAS_FUNC(dgemm)(
        (CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
        M, N, K, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

void amici_dgemv(
    BLASLayout layout, BLASTranspose TransA, int const M, int const N,
    double const alpha, double const* A, int const lda, double const* X,
    int const incX, double const beta, double* Y, int const incY
) {
    BLAS_FUNC(dgemv)(
        (CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A, lda, X,
        incX, beta, Y, incY
    );
}

void amici_daxpy(
    int const n, double const alpha, double const* x, int const incx, double* y,
    int const incy
) {
    BLAS_FUNC(daxpy)(n, alpha, x, incx, y, incy);
}

} // namespace amici
