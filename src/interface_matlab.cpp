/**
 * @file   interface_matlab.cpp
 * @brief  core routines for mex interface
 *
 * This file defines the function `mexFunction` which is executed upon calling
 * the `.mex` file from Matlab.
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

namespace amici {

    int dbl2int(const double x);

/**
 * @brief The mexRhsArguments enum takes care of the ordering of mex file
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
 * Translates AMICI_BLAS_TRANSPOSE values to CBLAS readable strings
 *
 * @param trans       flag indicating transposition and complex conjugation
 *
 * @return cblastrans CBlas readable CHAR indicating transposition and complex
 * conjugation
 */
char amici_blasCBlasTransToBlasTrans(BLASTranspose trans) {
    switch (trans) {
    case BLASTranspose::noTrans:
        return 'N';
    case BLASTranspose::trans:
        return 'T';
    case BLASTranspose::conjTrans:
        return 'C';
    }
    throw std::invalid_argument("Invalid argument to amici_blasCBlasTransToBlasTrans");
}

void amici_dgemm(BLASLayout layout, BLASTranspose TransA,
                 BLASTranspose TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    assert(layout == BLASLayout::colMajor);

    const ptrdiff_t M_ = M;
    const ptrdiff_t N_ = N;
    const ptrdiff_t K_ = K;
    const ptrdiff_t lda_ = lda;
    const ptrdiff_t ldb_ = ldb;
    const ptrdiff_t ldc_ = ldc;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);
    const char transB = amici_blasCBlasTransToBlasTrans(TransB);

    FORTRAN_WRAPPER(dgemm)(&transA, &transB,
                           &M_, &N_, &K_,
                           &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_);
}

void amici_dgemv(BLASLayout layout, BLASTranspose TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    assert(layout == BLASLayout::colMajor);

    const ptrdiff_t M_ = M;
    const ptrdiff_t N_ = N;
    const ptrdiff_t lda_ = lda;
    const ptrdiff_t incX_ = incX;
    const ptrdiff_t incY_ = incY;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);

    FORTRAN_WRAPPER(dgemv)(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
}

void amici_daxpy(int n, double alpha, const double *x, const int incx, double *y, int incy) {

    const ptrdiff_t n_ = n;
    const ptrdiff_t incx_ = incx;
    const ptrdiff_t incy_ = incy;

    FORTRAN_WRAPPER(daxpy)(&n_, &alpha, x, &incx_, y, &incy_);
}

/** conversion from mxArray to vector<realtype>
  * @param array Matlab array to create vector from
  * @param length Number of elements in array
  * @return std::vector with data from array
  */
std::vector<realtype> mxArrayToVector(const mxArray *array, int length) {
    return {mxGetPr(array), mxGetPr(array) + length};
}

std::unique_ptr<ExpData> expDataFromMatlabCall(const mxArray *prhs[],
                               Model const &model) {
    if (!mxGetPr(prhs[RHS_DATA]))
        return nullptr;

    auto edata = std::unique_ptr<ExpData>(new ExpData(model));

    // Y
    if (mxArray *dataY = mxGetProperty(prhs[RHS_DATA], 0, "Y")) {
        auto ny_my = static_cast<int>(mxGetN(dataY));
        if (ny_my != model.nytrue) {
            throw AmiException("Number of observables in data matrix (%i) does "
                               "not match model ny (%i)",
                               ny_my, model.nytrue);
        }
        auto nt_my = static_cast<int>(mxGetM(dataY));
        if (nt_my != model.nt()) {
            throw AmiException("Number of time-points in data matrix does (%i) "
                               "not match provided time vector (%i)",
                               nt_my, model.nt());
        }
        mxArray *dataYT;
        mexCallMATLAB(1, &dataYT, 1, &dataY, "transpose");
        auto observedData = mxArrayToVector(dataYT, ny_my * nt_my);
        edata->setObservedData(observedData);

    } else {
        throw AmiException("Field Y not specified as field in data struct!");
    }

    // Sigma Y
    if (mxArray *dataSigmaY = mxGetProperty(prhs[RHS_DATA], 0, "Sigma_Y")) {
        auto ny_sigmay = static_cast<int>(mxGetN(dataSigmaY));
        if (ny_sigmay != model.nytrue) {
            throw AmiException("Number of observables in data-sigma matrix (%i) "
                               "does not match model ny (%i)",
                               ny_sigmay, model.nytrue);
        }
        auto nt_sigmay = static_cast<int>(mxGetM(dataSigmaY));
        if (nt_sigmay != model.nt()) {
            throw AmiException("Number of time-points in data-sigma matrix (%i) "
                               "does not match provided time vector (%i)",
                               nt_sigmay, model.nt());
        }

        mxArray *dataSigmaYT;
        mexCallMATLAB(1, &dataSigmaYT, 1, &dataSigmaY, "transpose");
        auto observedDataSigma = mxArrayToVector(dataSigmaYT, ny_sigmay * nt_sigmay);
        edata->setObservedDataStdDev(observedDataSigma);
    } else {
        throw AmiException("Field Sigma_Y not specified as field in data struct!");
    }

    // Z
    if (mxArray *dataZ = mxGetProperty(prhs[RHS_DATA], 0, "Z")) {
        auto nz_mz = static_cast<int>(mxGetN(dataZ));
        if (nz_mz != model.nztrue) {
            throw AmiException("Number of events in event matrix (%i) does not "
                               "match provided nz (%i)",
                               nz_mz, model.nztrue);
        }
        auto ne_mz = static_cast<int>(mxGetM(dataZ));
        if (ne_mz != model.nMaxEvent()) {
            throw AmiException("Number of time-points in event matrix (%i) does "
                               "not match provided nmaxevent (%i)",
                               ne_mz, model.nMaxEvent());
        }
        mxArray *dataZT;
        mexCallMATLAB(1, &dataZT, 1, &dataZ, "transpose");
        auto observedEvents = mxArrayToVector(dataZT, nz_mz * ne_mz);
        edata->setObservedEvents(observedEvents);
    } else {
        throw AmiException("Field Z not specified as field in data struct!");
    }

    // Sigma Z
    if (mxArray *dataSigmaZ = mxGetProperty(prhs[RHS_DATA], 0, "Sigma_Z")) {
        auto nz_sigmaz = static_cast<int>(mxGetN(dataSigmaZ));
        if (nz_sigmaz != model.nztrue) {
            throw AmiException("Number of events in event-sigma matrix (%i) does "
                               "not match provided nz (%i)",
                               nz_sigmaz, model.nztrue);
        }
        auto ne_sigmaz = static_cast<int>(mxGetM(dataSigmaZ));
        if (ne_sigmaz != model.nMaxEvent()) {
            throw AmiException("Number of time-points in event-sigma matrix (%i) "
                               "does not match provided nmaxevent (%i)",
                               ne_sigmaz, model.nMaxEvent());
        }
        mxArray *dataSigmaZT;
        mexCallMATLAB(1, &dataSigmaZT, 1, &dataSigmaZ, "transpose");
        auto observedEventsSigma = mxArrayToVector(dataSigmaZT, nz_sigmaz * ne_sigmaz);
        edata->setObservedEventsStdDev(observedEventsSigma);
    } else {
        throw AmiException("Field Sigma_Z not specified as field in data struct!");

    }

    // preequilibration condition parameters
    if (mxArray *dataPreeq = mxGetProperty(prhs[RHS_DATA], 0, "conditionPreequilibration")) {
        int m = (int)mxGetM(dataPreeq);
        int n = (int)mxGetN(dataPreeq);
        if(m * n > 0) {
            if (m * n != model.nk() || (m != 1 && n != 1)) {
                throw AmiException("Number of preequilibration parameters (%dx%d) does "
                                   "not match model (%d)", m, n, model.nk());
            }
            edata->fixedParametersPreequilibration =
                    std::vector<realtype>(mxGetPr(dataPreeq), mxGetPr(dataPreeq) + m * n);
        }
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

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ss_atol")) {
            solver.setAbsoluteToleranceQuadratures(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ss_atol")));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ss_rtol")) {
            solver.setRelativeToleranceQuadratures(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ss_rtol")));
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
            solver.setSensitivityOrder(static_cast<SensitivityOrder>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ism")) {
            solver.setInternalSensitivityMethod(static_cast<InternalSensitivityMethod>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ism")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth")) {
            solver.setSensitivityMethod(static_cast<SensitivityMethod>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth_preeq")) {
            solver.setSensitivityMethodPreequilibration(
                static_cast<SensitivityMethod>(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "sensi_meth_preeq")))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "ordering")) {
            solver.setStateOrdering(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "ordering"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "stldet")) {
            solver.setStabilityLimitFlag(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "stldet"))));
        }

        if (mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_preeq")) {
            solver.setPreequilibration(dbl2int(mxGetScalar(mxGetProperty(prhs[RHS_OPTIONS], 0, "newton_preeq"))));
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
                model.setParameterScale(static_cast<ParameterScaling>(dbl2int(mxGetScalar(a))));
            } else if((mxGetM(a) == 1 && mxGetN(a) == model.np())
                      || (mxGetN(a) == 1 && mxGetM(a) == model.np())) {
                auto pscaleArray = static_cast<double *>(mxGetData(a));
                std::vector<ParameterScaling> pscale(model.np());
                for(int ip = 0; ip < model.np(); ++ip) {
                    pscale[ip] = static_cast<ParameterScaling>(dbl2int(pscaleArray[ip]));
                }
                model.setParameterScale(pscale);
            } else {
                throw AmiException("Provided pscale has invalid dimensions!");
            }
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
            if (mxGetM(x0) != model.nx_rdata) {
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
            if (mxGetM(x0) != model.nx_rdata) {
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
            if (mxGetM(sx0) != model.nx_rdata) {
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
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // use matlab error reporting
    amici::AmiciApplication amiciApp;
    amiciApp.warning = [](
            std::string const& identifier,
            std::string const& message){
        mexWarnMsgIdAndTxt(identifier.c_str(), message.c_str());
    };
    amiciApp.error = [](
            std::string const& identifier,
            std::string const& message){
        mexErrMsgIdAndTxt(identifier.c_str(), message.c_str());
    };

    if (nlhs != 1) {
        amiciApp.errorF("AMICI:mex:setup",
                        "Incorrect number of output arguments (must be 1)!");
    } else if(nrhs < amici::RHS_NUMARGS_REQUIRED) {
        amiciApp.errorF("AMICI:mex:setup",
                        "Incorrect number of input arguments (must be at least 7)!");
    };

    auto model = amici::generic_model::getModel();
    model->app = &amiciApp;

    auto solver = model->getSolver();
    solver->app = &amiciApp;
    setModelData(prhs, nrhs, *model);
    setSolverOptions(prhs, nrhs, *solver);

    std::unique_ptr<amici::ExpData> edata;
    if (nrhs > amici::RHS_DATA && mxGetPr(prhs[amici::RHS_DATA])) {
        try {
            edata = amici::expDataFromMatlabCall(prhs, *model);
        } catch (amici::AmiException const& ex) {
            amiciApp.errorF("AMICI:mex:setup","Failed to read experimental data:\n%s",ex.what());
        }
    } else if (solver->getSensitivityOrder() >= amici::SensitivityOrder::first &&
               solver->getSensitivityMethod() == amici::SensitivityMethod::adjoint) {
        amiciApp.errorF("AMICI:mex:setup","No data provided!");
    }

    /* ensures that plhs[0] is available */
    auto rdata = amici::runAmiciSimulation(*solver, edata.get(), *model);
    plhs[0] = getReturnDataMatlabFromAmiciCall(rdata.get());

}
