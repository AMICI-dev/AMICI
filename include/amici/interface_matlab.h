#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include "amici/amici.h"

#include <mex.h>
#include <memory>


namespace amici {

namespace generic_model {
extern std::unique_ptr<amici::Model> getModel();
} // namespace generic_model

class Model;
class ReturnDataMatlab;


/**
 * @brief setModelData sets data from the matlab call to the model object
 * @param prhs: pointer to the array of input arguments
 * @param nrhs: number of elements in prhs
 * @param model: model to update
 */
void setModelData(const mxArray *prhs[], int nrhs, Model& model);

/**
 * @brief setSolverOptions solver options from the matlab call to a solver
 * object
 * @param prhs: pointer to the array of input arguments
 * @param nrhs: number of elements in prhs
 * @param solver: solver to update
 */
void setSolverOptions(const mxArray *prhs[], int nrhs, Solver& solver);

/**
 * @brief setupReturnData initialises the return data struct
 * @param plhs user input
 * @param nlhs number of elements in plhs
 * @return rdata: return data struct
 */
ReturnDataMatlab *setupReturnData(mxArray *plhs[], int nlhs);


/*!
 * @brief expDataFromMatlabCall parses the experimental data from the matlab
 * call and writes it to an ExpData class object
 *
 * @param prhs pointer to the array of input arguments
 * @param model pointer to the model object, this is necessary to perform
 * dimension checks
 * @return edata pointer to experimental data object
 */
std::unique_ptr<ExpData> expDataFromMatlabCall(const mxArray *prhs[],
                                               const Model &model);

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
