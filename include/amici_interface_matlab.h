#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include <mex.h>
#include "include/amici.h"

class ReturnDataMatlab;

/**
 * @brief userDataFromMatlabCall extracts information from the matlab call and returns the corresponding UserData struct
 * @param[in] prhs: pointer to the array of input arguments @type mxArray
 * @return udata: struct containing all provided user data @type UserData
 */
UserData *userDataFromMatlabCall(const mxArray *prhs[]);

/**
 * setupReturnData initialises the return data struct
 * @param[in] plhs user input @type mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[out] pstatus pointer to the flag indicating the execution status @type double
 * @return rdata: return data struct @type ReturnData
 */
ReturnDataMatlab *setupReturnData(mxArray *plhs[], const UserData *udata, double *pstatus);

/**
 * expDataFromMatlabCall initialises the experimental data struct
 * @param[in] prhs user input @type *mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[out] status non-zero on failure, zero on success
 * @return edata: experimental data struct @type ExpData
 */
ExpData *expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata, int *status);


class ReturnDataMatlab : public ReturnData {

public:
    ReturnDataMatlab(const UserData *udata);
    ~ReturnDataMatlab() {}

    mxArray *mxsol;

protected:
    void initFields(const UserData *udata);

    virtual void initField1(double **fieldPointer, const char *fieldName, int dim);

    /**
     * @ brief initialise matrix and attach to the field
     * @ param FIELD name of the field to which the matrix will be attached
     * @ param D1 number of rows in the matrix
     * @ param D2 number of columns in the matrix
     */

    virtual void initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2);

    /**
     * @ brief initialise 3D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     */

    virtual void initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3);

    /**
     * @ brief initialise 4D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     * @ param D4 number of elements in the fourth dimension of the tensor
     */

    virtual void initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4);

};

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
