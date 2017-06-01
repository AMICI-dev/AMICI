#include "include/amici_interface_cpp.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <include/edata_accessors.h>
#include <include/tdata_accessors.h>

#include <cstring>

#define initField2(FIELD,D1,D2) \
FIELD ## data = new double[(D1) * (D2)]();
#define initField3(FIELD,D1,D2,D3) \
FIELD ## data = new double[(D1) * (D2) * (D3)]();

#define initField4(FIELD,D1,D2,D3,D4) \
FIELD ## data = new double[(D1) * (D2) * (D3) * (D4)]();


ReturnData *getSimulationResults(UserData *udata, const ExpData *edata) {
    double *originalParams = NULL;

    if(udata->pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * udata->np);
        memcpy(originalParams, udata->p, sizeof(double) * udata->np);
    }

    ReturnData *rdata = new ReturnData(udata);

    int status;
    runAmiciSimulation(udata, edata, rdata, &status);
    *rdata->status = status;

    if(originalParams) {
        memcpy(udata->p, originalParams, sizeof(double) * udata->np);
        free(originalParams);
    }

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
