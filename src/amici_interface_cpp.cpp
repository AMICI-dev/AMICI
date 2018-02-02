/**
 * @file   amici_interface_cpp.cpp
 * @brief  core routines for cpp interface
 *
 **/

#include "include/amici_interface_cpp.h"
#include "include/amici.h"
#include <include/amici_model.h>
#include <include/amici_exception.h>
#include <include/amici_solver.h>

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

#include <cstring>

namespace amici {

/*!
 * getSimulationResults is the core cpp interface function. It initializes the
 * model and return data and
 * then calls the core simulation routine.
 *
 * @param[in] model the model object, this is necessary to perform
 * dimension checks @type Model
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] solver solver object @type Solver
 * @return rdata pointer to return data object @type ReturnData
 */
ReturnData *getSimulationResults(Model &model, const ExpData *edata, Solver &solver) {

    ReturnData *rdata = new ReturnData(solver, &model);
    
    try {
        runAmiciSimulation(solver, edata, rdata, model);
        *rdata->status = AMICI_SUCCESS;
        rdata->applyChainRuleFactorToSimulationResults(&model);
    } catch (amici::IntegrationFailure& ex) {
        rdata->invalidate(ex.time);
        *(rdata->status) = ex.error_code;
    }

    return rdata;
}

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
void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
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
void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    cblas_dgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A,
                lda, X, incX, beta, Y, incY);
}

} // namespace amici
