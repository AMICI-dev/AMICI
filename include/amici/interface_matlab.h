#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include "amici/amici.h"

#include <mex.h>
#include <memory>

class Model;
extern std::unique_ptr<amici::Model> getModel();

namespace amici {    

class ReturnDataMatlab;


/**
 * @brief setModelData sets data from the matlab call to the model object
 * @param[in] prhs: pointer to the array of input arguments @type mxArray
 * @param[in] nrhs: number of elements in prhs
 * @param[in,out] model: model to update
 */
void setModelData(const mxArray *prhs[], int nrhs, Model& model);

/**
 * @brief setSolverOptions solver options from the matlab call to a solver object
 * @param[in] prhs: pointer to the array of input arguments @type mxArray
 * @param[in] nrhs: number of elements in prhs
 * @param[in,out] solver: solver to update
 */
void setSolverOptions(const mxArray *prhs[], int nrhs, Solver& solver);

/**
 * setupReturnData initialises the return data struct
 * @param[in] plhs user input @type mxArray
 * @param[in] nlhs number of elements in plhs @type mxArray
 * @return rdata: return data struct @type *ReturnData
 */
ReturnDataMatlab *setupReturnData(mxArray *plhs[], int nlhs);

/**
 * expDataFromMatlabCall initialises the experimental data struct
 * @param[in] prhs user input @type *mxArray
 * @return edata: experimental data struct @type *ExpData
 */
std::unique_ptr<ExpData> expDataFromMatlabCall(const mxArray *prhs[], const Model &model);

void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);

void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

} // namespace amici

#endif
