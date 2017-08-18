#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include <mex.h>
#include "include/amici.h"
#include "include/rdata.h"

class ReturnDataMatlab;

/**
 * @brief userDataFromMatlabCall extracts information from the matlab call and returns the corresponding UserData struct
 * @param[in] prhs: pointer to the array of input arguments @type mxArray
 * @return udata: struct containing all provided user data @type *UserData
 */
UserData *userDataFromMatlabCall(const mxArray *prhs[], int nrhs, Model *model);

/**
 * setupReturnData initialises the return data struct
 * @param[in] plhs user input @type mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[out] pstatus pointer to the flag indicating the execution status @type double
 * @return rdata: return data struct @type *ReturnData
 */
ReturnDataMatlab *setupReturnData(mxArray *plhs[], int nlhs, const UserData *udata);

/**
 * expDataFromMatlabCall initialises the experimental data struct
 * @param[in] prhs user input @type *mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[out] status non-zero on failure, zero on success
 * @return edata: experimental data struct @type *ExpData
 */
ExpData *expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata, Model *model);



void amici_dgemv(AMICI_BLAS_LAYOUT layout,
                 AMICI_BLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

#endif
