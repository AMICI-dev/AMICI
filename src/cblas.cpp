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

/*!
 * amici_dgemm provides an interface to the blas matrix matrix multiplication
 * routine dgemm. This routines computes
 * C = alpha*A*B + beta*C with A: [MxK] B:[KxN] C:[MxN]
 *
 * @param[in] layout    can be AMICI_BLAS_ColMajor or AMICI_BLAS_RowMajor.
 * @param[in] TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param[in] TransB    flag indicating whether B should be transposed before
 * multiplication
 * @param[in] M         number of rows in A/C
 * @param[in] N         number of columns in B/C
 * @param[in] K         number of rows in B, number of columns in A
 * @param[in] alpha     coefficient alpha
 * @param[in] A         matrix A
 * @param[in] lda       leading dimension of A (m or k)
 * @param[in] B         matrix B
 * @param[in] ldb       leading dimension of B (k or n)
 * @param[in] beta      coefficient beta
 * @param[in,out] C     matrix C
 * @param[in] ldc       leading dimension of C (m or n)
 */
void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    cblas_dgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA,
                (CBLAS_TRANSPOSE)TransB, M, N, K, alpha, A, lda, B, ldb, beta,
                C, ldc);
}

/*!
 * amici_dgemm provides an interface to the blas matrix vector multiplication
 * routine dgemv. This routines computes
 * y = alpha*A*x + beta*y with A: [MxK] B:[KxN] C:[MxN]
 *
 * @param[in] layout    can be AMICI_BLAS_ColMajor or AMICI_BLAS_RowMajor.
 * @param[in] TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param[in] M         number of rows in A
 * @param[in] N         number of columns in A
 * @param[in] alpha     coefficient alpha
 * @param[in] A         matrix A
 * @param[in] lda       leading dimension of A (m or n)
 * @param[in] X         vector X
 * @param[in] incX      increment for entries of X
 * @param[in] beta      coefficient beta
 * @param[in,out] Y     vector Y
 * @param[in] incY      increment for entries of Y
 */
void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    cblas_dgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A,
                lda, X, incX, beta, Y, incY);
}

/**
 * @brief Compute y = a*x + y
 * @param n number of elements in y
 * @param alpha scalar coefficient of x
 * @param x vector of length n*incx
 * @param incx x stride
 * @param y vector of length n*incy
 * @param incy y stride
 */
void amici_daxpy(int n, double alpha, const double *x, const int incx, double *y, int incy) {
    cblas_daxpy(n, alpha, x, incx, y, incy);
}

} // namespace amici
