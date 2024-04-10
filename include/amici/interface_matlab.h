#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include <amici/defines.h>

#include <memory>
#include <mex.h>

namespace amici {

class Model;
class Solver;
class ReturnDataMatlab;
class ExpData;

namespace generic_model {
extern std::unique_ptr<amici::Model> getModel();
} // namespace generic_model

/**
 * @brief setModelData sets data from the matlab call to the model object
 * @param prhs: pointer to the array of input arguments
 * @param nrhs: number of elements in prhs
 * @param model: model to update
 */
void setModelData(mxArray const* prhs[], int nrhs, Model& model);

/**
 * @brief setSolverOptions solver options from the matlab call to a solver
 * object
 * @param prhs: pointer to the array of input arguments
 * @param nrhs: number of elements in prhs
 * @param solver: solver to update
 */
void setSolverOptions(mxArray const* prhs[], int nrhs, Solver& solver);

/**
 * @brief setupReturnData initialises the return data struct
 * @param plhs user input
 * @param nlhs number of elements in plhs
 * @return rdata: return data struct
 */
ReturnDataMatlab* setupReturnData(mxArray* plhs[], int nlhs);

/*!
 * @brief expDataFromMatlabCall parses the experimental data from the matlab
 * call and writes it to an ExpData class object
 *
 * @param prhs pointer to the array of input arguments
 * @param model pointer to the model object, this is necessary to perform
 * dimension checks
 * @return edata pointer to experimental data object
 */
std::unique_ptr<ExpData>
expDataFromMatlabCall(mxArray const* prhs[], Model const& model);

void amici_dgemv(
    BLASLayout layout, BLASTranspose TransA, int const M, int const N,
    double const alpha, double const* A, int const lda, double const* X,
    int const incX, double const beta, double* Y, int const incY
);

void amici_dgemm(
    BLASLayout layout, BLASTranspose TransA, BLASTranspose TransB, int const M,
    int const N, int const K, double const alpha, double const* A,
    int const lda, double const* B, int const ldb, double const beta, double* C,
    int const ldc
);

} // namespace amici

#endif
