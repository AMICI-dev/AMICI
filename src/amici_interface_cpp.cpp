#include "include/amici_interface_cpp.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <include/edata_accessors.h>
#include <include/udata_accessors.h>
#include <include/rdata_accessors.h>
#include <include/tdata_accessors.h>

#include <cstring>

/**
 * @ brief initialise matrix and attach to the field
 * @ param FIELD name of the field to which the matrix will be attached
 * @ param D1 number of rows in the matrix
 * @ param D2 number of columns in the matrix
 */
#define initField2(FIELD,D1,D2) \
FIELD ## data = new double[(D1) * (D2)]();

/**
 * @ brief initialise 3D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 */
#define initField3(FIELD,D1,D2,D3) \
FIELD ## data = new double[(D1) * (D2) * (D3)]();

/**
 * @ brief initialise 4D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 * @ param D4 number of elements in the fourth dimension of the tensor
 */
#define initField4(FIELD,D1,D2,D3,D4) \
FIELD ## data = new double[(D1) * (D2) * (D3) * (D4)]();

ReturnData *initReturnData(const UserData *udata, int *pstatus) {
    ReturnData *rdata; /* returned rdata struct */

    /* Return rdata structure */
    rdata = new ReturnData();
    if (rdata == NULL)
        return(NULL);

    memset(rdata, 0, sizeof(*rdata));

    tsdata = new double[nt]();

    #include "include/amici_init_return_data_fields.h"

    return(rdata);
}


ReturnData *getSimulationResults(UserData *udata, const ExpData *edata, int *pstatus) {
    double *originalParams = NULL;

    if(udata->am_pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * np);
        memcpy(originalParams, p, sizeof(double) * np);
    }

    ReturnData *rdata = initReturnData(udata, pstatus);

    runAmiciSimulation(udata, edata, rdata, pstatus);

    if(originalParams) {
        memcpy(p, originalParams, sizeof(double) * np);
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
