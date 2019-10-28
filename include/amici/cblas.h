#ifndef AMICI_CBLAS_H
#define AMICI_CBLAS_H

#include "amici/defines.h"

namespace amici {

/**
 * amici_dgemm provides an interface to the CBlas matrix vector multiplication
 * routine dgemv. This routines computes
 * y = alpha*A*x + beta*y with A: [MxN] x:[Nx1] y:[Mx1]
 *
 * @param layout    always needs to be AMICI_BLAS_ColMajor.
 * @param TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param M         number of rows in A
 * @param N         number of columns in A
 * @param alpha     coefficient alpha
 * @param A         matrix A
 * @param lda       leading dimension of A (m or n)
 * @param X         vector X
 * @param incX      increment for entries of X
 * @param beta      coefficient beta
 * @param Y         vector Y
 * @param incY      increment for entries of Y
 */
void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 int M, int N, double alpha, const double *A,
                 int lda, const double *X, int incX,
                 double beta, double *Y, int incY);

/**
 * amici_dgemm provides an interface to the CBlas matrix matrix multiplication
 * routine dgemm. This routines computes
 * C = alpha*A*B + beta*C with A: [MxK] B:[KxN] C:[MxN]
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
 * @param lda       leading dimension of A (m or k)
 * @param B         matrix B
 * @param ldb       leading dimension of B (k or n)
 * @param beta      coefficient beta
 * @param C         matrix C
 * @param ldc       leading dimension of C (m or n)
 */
void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, int M, int N,
                 int K, double alpha, const double *A,
                 int lda, const double *B, int ldb,
                 double beta, double *C, int ldc);

/**
 * @brief Compute y = a*x + y
 * @param n         number of elements in y
 * @param alpha     scalar coefficient of x
 * @param x         vector of length n*incx
 * @param incx      x stride
 * @param y         vector of length n*incy
 * @param incy      y stride
 */
void amici_daxpy(int n, double alpha, const double *x, int incx, double *y, int incy);

} // namespace amici

#endif // AMICI_CBLAS_H
