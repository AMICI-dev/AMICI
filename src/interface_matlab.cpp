/**
 * @file   interface_matlab.cpp
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the
 * mex file from matlab
 */

#include "amici/interface_matlab.h"

#include "amici/model.h"
#include "amici/exception.h"
#include "amici/edata.h"
#include "amici/returndata_matlab.h"

#include <assert.h>
#include <blas.h>
#include <cstring>
#include <memory>

/** MS definition of PI and other constants */
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
/** define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

namespace amici {
    
    int dbl2int(const double x);

/**
 * @brief The mexFunctionArguments enum takes care of the ordering of mex file
 * arguments (indexing in prhs)
 */
enum mexRhsArguments {
    RHS_TIMEPOINTS,
    RHS_PARAMETERS,
    RHS_CONSTANTS,
    RHS_OPTIONS,
    RHS_PLIST,
    RHS_XSCALE_UNUSED,
    RHS_INITIALIZATION,
    RHS_DATA,
    RHS_NUMARGS_REQUIRED = RHS_DATA,
    RHS_NUMARGS
};


/*!
 * amici_blasCBlasTransToBlasTrans translates AMICI_BLAS_TRANSPOSE values to
 * CBlas readable strings
 *
 * @param trans       flag indicating transposition and complex conjugation
 *
 * @return cblastrans CBlas readable CHAR indicating transposition and complex
 * conjugation
 */
char amici_blasCBlasTransToBlasTrans(AMICI_BLAS_TRANSPOSE trans) {
    switch (trans) {
    case AMICI_BLAS_NoTrans:
        return 'N';
    case AMICI_BLAS_Trans:
        return 'T';
    case AMICI_BLAS_ConjTrans:
        return 'C';
    }
}

/*!
 * amici_dgemm provides an interface to the CBlas matrix matrix multiplication
 * routine dgemm. This routines computes
 * C = alpha*A*B + beta*C with A: [MxK] B:[KxN] C:[MxN]
 *
 * @param layout    memory layout.
 * @param TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param TransB    flag indicating whether B should be transposed before
 * multiplication
 * @param M         number of rows in A/C
 * @param N         number of columns in B/C
 * @param K         number of rows in B, number of columns in A
 * @param alpha     coefficient alpha
 * @param A         matrix A
 * @param lda       leading dimension of A (m or k)
 * @param B         matrix B
 * @param ldb       leading dimension of B (k or n)
 * @param beta      coefficient beta
 * @param[in,out] C     matrix C
 * @param ldc       leading dimension of C (m or n)
 * @return void
 */
void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    //assert(layout == AMICI_BLAS_RowMajor);

    const ptrdiff_t M_ = M;
    const ptrdiff_t N_ = N;
    const ptrdiff_t K_ = K;
    const ptrdiff_t lda_ = lda;
    const ptrdiff_t ldb_ = ldb;
    const ptrdiff_t ldc_ = ldc;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);
    const char transB = amici_blasCBlasTransToBlasTrans(TransB);

#if defined(_WIN32)
    dgemm(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C,
          &ldc_);
#else
    dgemm_(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta,
           C, &ldc_);
#endif
}

/*!
 * amici_dgemm provides an interface to the CBlas matrix vector multiplication
 * routine dgemv. This routines computes
 * y = alpha*A*x + beta*y with A: [MxN] x:[Nx1] y:[Mx1]
 *
 * @param layout    always needs to be AMICI_BLAS_ColMajor.
 * @param TransA    flag indicating whether A should be transposed before
 * multiplication
 * @param M         number of rows in A
 * @param N         number of columns in A
 * @param alpha     coefficient alpha
 * @param A         matrix A
 * @param lda       leading dimension of A (m or n)
 * @param X         vector X
 * @param incX      increment for entries of X
 * @param beta      coefficient beta
 * @param[in,out] Y     vector Y
 * @param incY      increment for entries of Y
 * @return void
 */
void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    assert(layout == AMICI_BLAS_RowMajor);

    const ptrdiff_t M_ = M;
    const ptrdiff_t N_ = N;
    const ptrdiff_t lda_ = lda;
    const ptrdiff_t incX_ = incX;
    const ptrdiff_t incY_ = incY;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);

#if defined(_WIN32)
    dgemv(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#else
    dgemv_(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#endif
}

/*!
 * expDataFromMatlabCall parses the experimental data from the matlab call and
 * writes it to an ExpData class object
 *
 * @param prhs pointer to the array of input arguments
 * @param model pointer to the model object, this is necessary to perform
 * dimension checks
 * @return edata pointer to experimental data object
 */
ExpData *expDataFromMatlabCall(const mxArray *prhs[],
                               Model const &model) {
    if (!mxGetPr(prhs[RHS_DATA]))
        return nullptr;

    int nt_my = 0, ny_my = 0, nt_sigmay = 0,
        ny_sigmay = 0; /* integers with problem dimensionality */
    int ne_mz = 0, nz_mz = 0, ne_sigmaz = 0,
        nz_sigmaz = 0; /* integers with problem dimensionality */
    ExpData *edata = new ExpData(model);
    if (edata->my.empty() && edata->mz.empty()) {
        // if my and mz are both empty, no (or empty) data was provided
        // in that case we simply return a nullptr for easier checking.
        delete(edata);
        return nullptr;
    }
    
    // Y
    if (mxArray *dataY = mxGetProperty(prhs[RHS_DATA], 0, "Y")) {
        ny_my = (int)mxGetN(dataY);
        if (ny_my != model.nytrue) {
            throw AmiException("Number of observables in data matrix (%i) does "
                               "not match model ny (%i)",
                               ny_my, model.nytrue);
        }
        nt_my = (int)mxGetM(dataY);
        if (nt_my != model.nt()) {
            throw AmiException("Number of time-points in data matrix does (%i) "
                               "not match provided time vector (%i)",
                               nt_my, model.nt());
        }
        mxArray *dataYT;
        mexCallMATLAB(1, &dataYT, 1, &dataY, "transpose");
        edata->setObservedData(mxGetPr(dataYT));
        
    } else {
        throw AmiException("Field Y not specified as field in data struct!");
    }
    
    // Sigma Y
    if (mxArray *dataSigmaY = mxGetProperty(prhs[RHS_DATA], 0, "Sigma_Y")) {
        ny_sigmay = (int)mxGetN(dataSigmaY);
        if (ny_sigmay != model.nytrue) {
            throw AmiException("Number of observables in data-sigma matrix (%i) "
                               "does not match model ny (%i)",
                               ny_sigmay, model.nytrue);
        }
        nt_sigmay = (int)mxGetM(dataSigmaY);
        if (nt_sigmay != model.nt()) {
            throw AmiException("Number of time-points in data-sigma matrix (%i) "
                               "does not match provided time vector (%i)",
                               nt_sigmay, model.nt());
        }
        
        mxArray *dataSigmaYT;
        mexCallMATLAB(1, &dataSigmaYT, 1, &dataSigmaY, "transpose");
        edata->setObservedDataStdDev(mxGetPr(dataSigmaYT));
    } else {
        throw AmiException("Field Sigma_Y not specified as field in data struct!");
    }
    
    // Z
    if (mxArray *dataZ = mxGetProperty(prhs[RHS_DATA], 0, "Z")) {
        nz_mz = (int)mxGetN(dataZ);
        if (nz_mz != model.nztrue) {
            throw AmiException("Number of events in event matrix (%i) does not "
                               "match provided nz (%i)",
                               nz_mz, model.nztrue);
        }
        ne_mz = (int)mxGetM(dataZ);
        if (ne_mz != model.nMaxEvent()) {
            throw AmiException("Number of time-points in event matrix (%i) does "
                               "not match provided nmaxevent (%i)",
                               ne_mz, model.nMaxEvent());
        }
        mxArray *dataZT;
        mexCallMATLAB(1, &dataZT, 1, &dataZ, "transpose");
        edata->setObservedEvents(mxGetPr(dataZT));
    } else {
        throw AmiException("Field Z not specified as field in data struct!");
    }
    
    // Sigma Z
    if (mxArray *dataSigmaZ = mxGetProperty(prhs[RHS_DATA], 0, "Sigma_Z")) {
        nz_sigmaz = (int)mxGetN(dataSigmaZ);
        if (nz_sigmaz != model.nztrue) {
            throw AmiException("Number of events in event-sigma matrix (%i) does "
                               "not match provided nz (%i)",
                               nz_sigmaz, model.nztrue);
        }
        ne_sigmaz = (int)mxGetM(dataSigmaZ);
        if (ne_sigmaz != model.nMaxEvent()) {
            throw AmiException("Number of time-points in event-sigma matrix (%i) "
                               "does not match provided nmaxevent (%i)",
                               ne_sigmaz, model.nMaxEvent());
        }
        mxArray *dataSigmaZT;
        mexCallMATLAB(1, &dataSigmaZT, 1, &dataSigmaZ, "transpose");
        edata->setObservedEventsStdDev(mxGetPr(dataSigmaZT));
    } else {
        throw AmiException("Field Sigma_Z not specified as field in data struct!");
        
    }
    return edata;
}

/** conversion from double to int with checking for loss of data
  *  @param x input
  *  @return int_x casted value
  */
int dbl2int(const double x){
    if((std::round(x)-x) != 0.0)
        throw AmiException("Invalid non-integer value for integer option");
    return(static_cast<int>(x));
}

void setSolverOptions(const mxArray *prhs[], int nrhs, Solver &solver)
{
    if (mxGetPr(prhs[RHS_OPTIONS])) {
        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "atol")) {
            solver.setAbsoluteTolerance(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "atol")));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "rtol")) {
            solver.setRelativeTolerance(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "rtol")));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "quad_atol")) {
            solver.setAbsoluteToleranceQuadratures(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "quad_atol")));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "quad_rtol")) {
            solver.setRelativeToleranceQuadratures(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "quad_rtol")));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "maxsteps")) {
            solver.setMaxSteps(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "maxsteps"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "maxstepsB")) {
            solver.setMaxStepsBackwardProblem(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "maxstepsB"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "lmm")) {
            solver.setLinearMultistepMethod(static_cast<LinearMultistepMethod>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "lmm")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "iter")) {
            solver.setNonlinearSolverIteration(static_cast<NonlinearSolverIteration>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "iter")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "interpType")) {
            solver.setInterpolationType(static_cast<InterpolationType>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "interpType")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "linsol")) {
            solver.setLinearSolver(static_cast<LinearSolver>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "linsol")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi")) {
            solver.setSensitivityOrder(static_cast<AMICI_sensi_order>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ism")) {
            solver.setInternalSensitivityMethod(static_cast<InternalSensitivityMethod>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ism")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth")) {
            solver.setSensitivityMethod(static_cast<AMICI_sensi_meth>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ordering")) {
            solver.setStateOrdering(static_cast<StateOrdering>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ordering")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "stldet")) {
            solver.setStabilityLimitFlag(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "stldet"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_preeq")) {
            solver.setNewtonPreequilibration(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_preeq"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_precon")) {
            solver.setNewtonPreconditioner(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_precon"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_maxsteps")) {
            solver.setNewtonMaxSteps(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_maxsteps"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_maxlinsteps")) {
            solver.setNewtonMaxLinearSteps(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_maxlinsteps"))));
        }
    }
}

void setModelData(const mxArray *prhs[], int nrhs, Model &model)
{
    if (mxGetPr(prhs[RHS_OPTIONS])) {
        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "nmaxevent")) {
            model.setNMaxEvent(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "nmaxevent"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "tstart")) {
            model.setT0(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "tstart")));
        }

        if (mxArray *a = mxGetProperty(prhs[RHS_OPTIONS], 0, "pscale")) {
            if(mxGetM(a) == 1 && mxGetN(a) == 1) {
                model.setParameterScale(static_cast<AMICI_parameter_scaling>(dbl2int(mxGetScalar(a))));
            } else if((mxGetM(a) == 1 && mxGetN(a) == model.np())
                      || (mxGetN(a) == 1 && mxGetM(a) == model.np())) {
                auto pscaleArray = static_cast<double *>(mxGetData(a));
                std::vector<AMICI_parameter_scaling> pscale(model.np());
                for(int ip = 0; ip < model.np(); ++ip) {
                    pscale[ip] = static_cast<AMICI_parameter_scaling>(dbl2int(pscaleArray[ip]));
                }
                model.setParameterScale(pscale);
            } else {
                throw AmiException("Provided pscale has invalid dimensions!");
            }
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "qpositivex")) {
            mxArray *a = mxGetProperty(prhs[RHS_OPTIONS], 0, "qpositivex");
            int len = (int)mxGetM(a) * mxGetN(a);
            if(mxGetM(a) != 1 && mxGetN(a) != 1)
                throw AmiException("Provided qpositivex has invalid dimensions!");
            model.setPositivityFlag(std::vector<int>((double *)mxGetData(a),(double *)mxGetData(a)+len));
        }
    }

    if (prhs[RHS_TIMEPOINTS] &&
                mxGetM(prhs[RHS_TIMEPOINTS]) * mxGetN(prhs[RHS_TIMEPOINTS]) > 0) {
            model.setTimepoints(std::vector<double>(
                                    mxGetPr(prhs[RHS_TIMEPOINTS]),
                                    mxGetPr(prhs[RHS_TIMEPOINTS])
                                    + (int)mxGetM(prhs[RHS_TIMEPOINTS]) * mxGetN(prhs[RHS_TIMEPOINTS])));

    }

    if (model.np() > 0) {
        if (mxGetPr(prhs[RHS_PARAMETERS])) {
            if (mxGetM(prhs[RHS_PARAMETERS]) * mxGetN(prhs[RHS_PARAMETERS]) ==
                    model.np()) {
                model.setParameters(std::vector<double>(mxGetPr(prhs[RHS_PARAMETERS]),
                                                        mxGetPr(prhs[RHS_PARAMETERS])
                                                        + mxGetM(prhs[RHS_PARAMETERS]) * mxGetN(prhs[RHS_PARAMETERS])));
            }
        }
    }

    if (model.nk() > 0) {
        if (mxGetPr(prhs[RHS_CONSTANTS])) {
            if (mxGetM(prhs[RHS_CONSTANTS]) * mxGetN(prhs[RHS_CONSTANTS]) ==
                    model.nk()) {
                model.setFixedParameters(std::vector<double>(mxGetPr(prhs[RHS_CONSTANTS]),
                                                             mxGetPr(prhs[RHS_CONSTANTS])
                                                             + mxGetM(prhs[RHS_CONSTANTS]) * mxGetN(prhs[RHS_CONSTANTS])));
            }
        }
    }
    if (mxGetPr(prhs[RHS_PLIST])) {
        model.setParameterList(std::vector<int>(mxGetPr(prhs[RHS_PLIST]),
                                                mxGetPr(prhs[RHS_PLIST])
                                                + mxGetM(prhs[RHS_PLIST]) * mxGetN(prhs[RHS_PLIST])));
    } else {
        model.requireSensitivitiesForAllParameters();
    }

    /* Check, if initial states and sensitivities are passed by user or must be
             * calculated */
    if (mxGetPr(prhs[RHS_INITIALIZATION])) {
        mxArray *x0 = mxGetField(prhs[RHS_INITIALIZATION], 0, "x0");
        if (x0 && (mxGetM(x0) * mxGetN(x0)) > 0) {
            /* check dimensions */
            if (mxGetN(x0) != 1) {
                throw AmiException("Number of rows in x0 field must be equal to 1!");
            }
            if (mxGetM(x0) != model.nx) {
                throw AmiException("Number of columns in x0 field "
                                   "does not agree with number of "
                                   "model states!");
            }
        }
    }

    /* Check, if initial states and sensitivities are passed by user or must be
             * calculated */
    if (mxGetPr(prhs[RHS_INITIALIZATION])) {
        mxArray *x0 = mxGetField(prhs[RHS_INITIALIZATION], 0, "x0");
        if (x0 && (mxGetM(x0) * mxGetN(x0)) > 0) {
            /* check dimensions */
            if (mxGetN(x0) != 1) {
                throw AmiException("Number of rows in x0 field must be equal to 1!");
            }
            if (mxGetM(x0) != model.nx) {
                throw AmiException("Number of columns in x0 field "
                                   "does not agree with number of "
                                   "model states!");
            }

            model.setInitialStates(std::vector<double>(mxGetPr(x0),
                                                       mxGetPr(x0) + mxGetM(x0) * mxGetN(x0)));
        }

        mxArray *sx0 = mxGetField(prhs[RHS_INITIALIZATION], 0, "sx0");
        if (sx0 && (mxGetM(sx0) * mxGetN(sx0)) > 0) {
            /* check dimensions */
            if (mxGetN(sx0) != model.nplist()) {
                throw AmiException("Number of rows in sx0 field "
                                   "does not agree with number of "
                                   "model parameters!");
            }
            if (mxGetM(sx0) != model.nx) {
                throw AmiException("Number of columns in sx0 "
                                   "field does not agree with "
                                   "number of model states!");
            }
            model.setInitialStateSensitivities(std::vector<double>(mxGetPr(sx0),
                                                                   mxGetPr(sx0) + mxGetM(sx0) * mxGetN(sx0)));
        }
    }
}

} // namespace amici


/*!
 * mexFunction is the main interface function for the MATLAB interface. It reads
 * in input data (udata and edata) and
 * creates output data compound (rdata) and then calls the AMICI simulation
 * routine to carry out numerical integration.
 *
 * @param nlhs number of output arguments of the matlab call
 * @param plhs pointer to the array of output arguments
 * @param nrhs number of input arguments of the matlab call
 * @param prhs pointer to the array of input arguments
 * @return void
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // use matlab error reporting
    amici::warnMsgIdAndTxt = &mexWarnMsgIdAndTxt;
    amici::errMsgIdAndTxt = &mexErrMsgIdAndTxt;
    
    if (nlhs != 1) {
        amici::errMsgIdAndTxt("AMICI:mex:setup","Incorrect number of output arguments (must be 1)!");
    } else if(nrhs < amici::RHS_NUMARGS_REQUIRED) {
        amici::errMsgIdAndTxt("AMICI:mex:setup", "Incorrect number of input arguments (must be at least 7)!");
    };

    auto model = getModel();
    auto solver = model->getSolver();
    setModelData(prhs, nrhs, *model);
    setSolverOptions(prhs, nrhs, *solver);

    std::unique_ptr<amici::ExpData> edata;
    if (nrhs > amici::RHS_DATA && mxGetPr(prhs[amici::RHS_DATA])) {
        try {
            edata.reset(amici::expDataFromMatlabCall(prhs, *model));
        } catch (amici::AmiException& ex) {
            amici::errMsgIdAndTxt("AMICI:mex:setup","Failed to read experimental data:\n%s",ex.what());
        }
    } else if (solver->getSensitivityOrder() >= amici::AMICI_SENSI_ORDER_FIRST &&
               solver->getSensitivityMethod() == amici::AMICI_SENSI_ASA) {
        amici::errMsgIdAndTxt("AMICI:mex:setup","No data provided!");
    }
    
    /* ensures that plhs[0] is available */
    auto rdata = amici::runAmiciSimulation(*solver, edata.get(), *model);
    plhs[0] = getReturnDataMatlabFromAmiciCall(rdata.get());
    
}


