#ifndef AMICI_CBLAS_H
#define AMICI_CBLAS_H

#include "amici/defines.h"

namespace amici {

/**
 * @brief CBLAS matrix vector multiplication (dgemv).
 *
 * Computes \f$ y = alpha*A*x + beta*y \f$ with A: [MxN] x:[Nx1] y:[Mx1]
 *
 * @param layout    Matrix layout, column major or row major.
 * @param TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param M         number of rows in A
 * @param N         number of columns in A
 * @param alpha     coefficient alpha
 * @param A         matrix A
 * @param lda       leading dimension / stride of A (>=N if row-major,
 * >=M if col-major)
 * @param X         vector X
 * @param incX      increment for entries of X
 * @param beta      coefficient beta
 * @param Y         vector Y
 * @param incY      increment for entries of Y
 */
void amici_dgemv(
    BLASLayout layout, BLASTranspose TransA, int M, int N, double alpha,
    double const* A, int lda, double const* X, int incX, double beta, double* Y,
    int incY
);

/**
 * @brief CBLAS matrix matrix multiplication (dgemm)
 *
 * This routines computes \f$ C = alpha*A*B + beta*C \f$
 * with A: [MxK] B:[KxN] C:[MxN]
 *
 * @param layout    memory layout.
 * @param TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param TransB    flag indicating whether B should be transposed before
 * multiplication
 * @param M         number of rows in A/C
 * @param N         number of columns in B/C
 * @param K         number of rows in B, number of columns in A
 * @param alpha     coefficient alpha
 * @param A         matrix A
 * @param lda       leading dimension of A (>=M or >=K)
 * @param B         matrix B
 * @param ldb       leading dimension of B (>=K or >=N)
 * @param beta      coefficient beta
 * @param C         matrix C
 * @param ldc       leading dimension of C (>=M or >= N)
 */
void amici_dgemm(
    BLASLayout layout, BLASTranspose TransA, BLASTranspose TransB, int M, int N,
    int K, double alpha, double const* A, int lda, double const* B, int ldb,
    double beta, double* C, int ldc
);

/**
 * @brief Compute y = a*x + y
 * @param n         number of elements in y
 * @param alpha     scalar coefficient of x
 * @param x         vector of length n*incx
 * @param incx      x stride
 * @param y         vector of length n*incy
 * @param incy      y stride
 */
void amici_daxpy(
    int n, double alpha, double const* x, int incx, double* y, int incy
);

} // namespace amici

#endif // AMICI_CBLAS_H
