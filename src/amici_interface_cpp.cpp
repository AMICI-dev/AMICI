#include "include/amici_interface_cpp.h"
#include <include/amici_model.h>
#include <include/amici_model_functions.h>
#include <include/amici_solver.h>
#include "include/amici.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif defined(AMICI_BLAS_MKL)
#include <mkl.h>
#else
#include <cblas.h>
#endif

#include <cstring>

ReturnData *getSimulationResults(UserData *udata, const ExpData *edata) {

    Model *model = getModel();
    Solver *solver = getSolver();

    ReturnData *rdata = new ReturnData(udata, model);

    int status = runAmiciSimulation(udata, edata, rdata,  model, solver);
    *rdata->status = status;

    delete model;
    delete solver;

    return rdata;
}

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA, AMICI_BLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
    cblas_dgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta, double *Y, const int incY)
{
    cblas_dgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
