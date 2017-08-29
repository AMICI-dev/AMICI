/**
 * @file   amici_interface_matlab.cpp
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the
 * mex file from matlab
 */

#include "include/amici_interface_matlab.h"
#include "include/amici_model.h"
#include "include/edata.h"
#include "include/returndata_matlab.h"
#include "include/udata.h"

#include <assert.h>
#include <blas.h>
#include <cstring>

/** MS definition of PI and other constants */
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI 
/** define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

/**
 * @ brief extract information from a property of a matlab class (scalar)
 * @ param OPTION name of the property
 * @ param TYPE class to which the information should be cast
 */
#define readOptionScalar(OPTION, TYPE)                                         \
    if (mxGetProperty(prhs[3], 0, #OPTION)) {                                  \
        udata->OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3], 0, #OPTION)); \
    } else {                                                                   \
        warnMsgIdAndTxt("AMICI:mex:OPTION",                                    \
                        "Provided options do not have field " #OPTION "!");    \
        goto freturn;                                                          \
    }

/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION)                                                 \
    if (mxGetProperty(prhs[3], 0, #OPTION)) {                                  \
        mxArray *a = mxGetProperty(prhs[3], 0, #OPTION);                       \
        int len = (int)mxGetM(a) * mxGetN(a);                                  \
        udata->OPTION = new double[len];                                       \
        memcpy(udata->OPTION, mxGetData(a), sizeof(double) * len);             \
    } else {                                                                   \
        warnMsgIdAndTxt("AMICI:mex:OPTION",                                    \
                        "Provided options do not have field " #OPTION "!");    \
        goto freturn;                                                          \
    }

/*!
 * mexFunction is the main interface function for the MATLAB interface. It reads in input data (udata and edata) and 
 * creates output data compound (rdata) and then calls the AMICI simulation routine to carry out numerical integration.
 *
 * @param[in] nlhs number of output arguments of the matlab call @type int
 * @param[out] plhs pointer to the array of output arguments @type mxArray
 * @param[in] nrhs number of input arguments of the matlab call @type int
 * @param[in] prhs pointer to the array of input arguments @type mxArray
 * @return void
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // use matlab error reporting
    warnMsgIdAndTxt = &mexWarnMsgIdAndTxt;
    errMsgIdAndTxt = &mexErrMsgIdAndTxt;

    /* ensures that plhs[0] is available */
    if (nlhs != 1) {
        errMsgIdAndTxt("AMICI:mex",
                       "Incorrect number of output arguments (must be 1)!");
        return;
    };

    Model *model = getModel();

    UserData *udata = userDataFromMatlabCall(prhs, nrhs, model);

    ReturnDataMatlab *rdata = new ReturnDataMatlab(udata, model);
    plhs[0] = rdata->mxsol;
    if (*(rdata->status) != AMICI_SUCCESS)
        return;

    ExpData *edata = NULL;
    if (nrhs >= 8) {
        edata = expDataFromMatlabCall(prhs, udata, model);
    }

    if (udata) {
        if (model->nx > 0) {
            *(rdata->status) =
                (double)runAmiciSimulation(udata, edata, rdata, model);
        }
    }

    delete model;

    if (udata)
        delete udata;
    if (rdata)
        delete rdata;
    if (edata)
        delete edata;
}

/*!
 * userDataFromMatlabCall parses the input from the matlab call and writes it to an UserData class object
 *
 * @param[in] nrhs number of input arguments of the matlab call @type int
 * @param[in] prhs pointer to the array of input arguments @type mxArray
 * @param[in] model pointer to the model object, this is necessary to perform dimension checks @type Model
 * @return udata pointer to user data object @type UserData
 */
UserData *userDataFromMatlabCall(const mxArray *prhs[], int nrhs,
                                 Model *model) {
    if (nrhs < 8) {
        errMsgIdAndTxt(
            "AMICI:mex",
            "Incorrect number of input arguments (must be at least 7)!");
        return NULL;
    };

    UserData *udata = new UserData();

    /* time */
    if (prhs[0] || mxGetM(prhs[0]) * mxGetN(prhs[0]) == 0) {
        udata->nt = (int)mxGetM(prhs[0]) * mxGetN(prhs[0]);
        udata->ts = new double[udata->nt];
        memcpy(udata->ts, mxGetPr(prhs[0]), sizeof(double) * udata->nt);
    } else {
        errMsgIdAndTxt("AMICI:mex:tout", "No time vector provided!");
        goto freturn;
    }

    /*
     `if(udata->nX > 0)` checks here are necessary as respective [1x0] MATLAB
     inputs
     will procude NULL pointers when calling mxGetPr. checking nrhs ensures that
     prhs[irhs]
     actually exist, otherwise mxGetPr will segfault. (checking phrs[irhs] does
     not help)
     */

    /* parameters */

    if (model->np > 0) {
        if (mxGetPr(prhs[1])) {
            if (mxGetM(prhs[1]) * mxGetN(prhs[1]) == model->np) {
                udata->p = new double[model->np];
                memcpy(udata->p, mxGetPr(prhs[1]), sizeof(double) * model->np);
            } else {
                errMsgIdAndTxt(
                    "AMICI:mex:theta",
                    "Provided parameter vector has incorrect length!");
                goto freturn;
            }
        } else {
            errMsgIdAndTxt("AMICI:mex:theta", "No parameter vector provided!");
            goto freturn;
        }
    }

    /* constants */
    if (model->nk > 0) {
        if (mxGetPr(prhs[2])) {
            if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == model->nk) {
                udata->k = new double[model->nk];
                memcpy(udata->k, mxGetPr(prhs[2]), sizeof(double) * model->nk);
            } else {
                errMsgIdAndTxt(
                    "AMICI:mex:kappa",
                    "Provided constant vector has incorrect length!");
                goto freturn;
            }
        } else {
            errMsgIdAndTxt("AMICI:mex:kappa", "No constant vector provided!");
            goto freturn;
        }
    }

    /* options */
    if (mxGetPr(prhs[3])) {
        readOptionScalar(nmaxevent, int) readOptionScalar(tstart, double)
            readOptionScalar(atol, double) readOptionScalar(rtol, double)
                readOptionScalar(maxsteps, int) readOptionScalar(lmm, int)
                    readOptionScalar(iter, int) readOptionScalar(
                        interpType, int) readOptionScalar(linsol, int)
                        readOptionScalar(stldet, booleantype) readOptionData(
                            qpositivex) readOptionScalar(sensi,
                                                         AMICI_sensi_order)
                            readOptionScalar(pscale, AMICI_parameter_scaling)
                                readOptionScalar(ism, int) readOptionScalar(
                                    sensi_meth, AMICI_sensi_meth)
                                    readOptionScalar(ordering, int)
                                        readOptionScalar(newton_preeq, int)
                                            readOptionScalar(newton_precon, int)
                                                readOptionScalar(
                                                    newton_maxsteps, int)
                                                    readOptionScalar(
                                                        newton_maxlinsteps, int)
    } else {
        errMsgIdAndTxt("AMICI:mex:options", "No options provided!");
        delete udata;
        goto freturn;
    }

    /* plist */
    if (mxGetPr(prhs[4])) {
        udata->nplist = (int)mxGetM(prhs[4]) * mxGetN(prhs[4]);
        udata->plist = new int[udata->nplist]();
        realtype *plistdata = mxGetPr(prhs[4]);

        for (int ip = 0; ip < udata->nplist; ip++) {
            udata->plist[ip] = (int)plistdata[ip];
        }
    } else {
        if (udata->sensi != AMICI_SENSI_ORDER_NONE) {
            errMsgIdAndTxt("AMICI:mex:plist", "No parameter list provided!");
            goto freturn;
        } else {
            udata->nplist = 0;
        }
    }

    /* pbar */
    if (udata->nplist > 0) {
        if (mxGetPr(prhs[5])) {
            if (mxGetM(prhs[5]) * mxGetN(prhs[5]) == udata->nplist) {
                udata->pbar = new double[udata->nplist]();
                memcpy(udata->pbar, mxGetPr(prhs[5]),
                       sizeof(double) * udata->nplist);
            } else {
                errMsgIdAndTxt(
                    "AMICI:mex:pbar",
                    "Provided parameter scales have incorrect length!");
                goto freturn;
            }
        } else {
            errMsgIdAndTxt("AMICI:mex:pbar", "No parameter scales provided!");
            goto freturn;
        }
    }

    /* xscale */
    /* this check previously always failed (xscale was empty all along).
    Commented out until implemented
    if (mxGetPr(prhs[6])) {
        if( mxGetM(prhs[6]) * mxGetN(prhs[6]) == udata->nx) {
            udata->xbar = new double[udata->nx]();
            memcpy(udata->xbar, mxGetPr(prhs[6]), sizeof(double) * udata->nx);
        } else {
            errMsgIdAndTxt("AMICI:mex:xbar","Provided state scales have
    incorrect length!");
            goto freturn;
        }
    } else {
        errMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
        goto freturn;
    }*/

    /* Check, if initial states and sensitivities are passed by user or must be
     * calculated */
    if (mxGetPr(prhs[7])) {
        mxArray *x0 = mxGetField(prhs[7], 0, "x0");
        if (x0 && (mxGetM(x0) * mxGetN(x0)) > 0) {
            /* check dimensions */
            if (mxGetN(x0) != 1) {
                errMsgIdAndTxt(
                    "AMICI:mex:x0",
                    "Number of rows in x0 field must be equal to 1!");
                goto freturn;
            }
            if (mxGetM(x0) != model->nx) {
                errMsgIdAndTxt("AMICI:mex:x0", "Number of columns in x0 field "
                                               "does not agree with number of "
                                               "model states!");
                goto freturn;
            }

            if ((mxGetM(x0) * mxGetN(x0)) > 0) {
                udata->x0data = new double[model->nx];
                memcpy(udata->x0data, mxGetPr(x0), sizeof(double) * model->nx);
            }
        }

        mxArray *sx0 = mxGetField(prhs[7], 0, "sx0");
        if (sx0 && (mxGetM(sx0) * mxGetN(sx0)) > 0) {
            /* check dimensions */
            if (mxGetN(sx0) != udata->nplist) {
                errMsgIdAndTxt("AMICI:mex:sx0", "Number of rows in sx0 field "
                                                "does not agree with number of "
                                                "model parameters!");
                goto freturn;
            }
            if (mxGetM(sx0) != model->nx) {
                errMsgIdAndTxt("AMICI:mex:sx0", "Number of columns in sx0 "
                                                "field does not agree with "
                                                "number of model states!");
                goto freturn;
            }

            udata->sx0data = new double[model->nx * udata->nplist];
            memcpy(udata->sx0data, mxGetPr(sx0),
                   sizeof(double) * model->nx * udata->nplist);
        }
    }
    return udata;

freturn:
    delete udata;
    return NULL;
}


/*!
 * amici_blasCBlasTransToBlasTrans translates AMICI_BLAS_TRANSPOSE values to CBlas readable strings
 *
 * @param[in] trans       flag indicating transposition and complex conjugation @type AMICI_BLAS_TRANSPOSE
 * @return cblastrans CBlas readable CHAR indicating transposition and complex conjugation @type char
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
 * amici_dgemm provides an interface to the CBlas matrix matrix multiplication routine dgemm. This routines computes
 * C = alpha*A*B + beta*C with A: [MxK] B:[KxN] C:[MxN]
 *
 * @param[in] layout    always needs to be AMICI_BLAS_ColMajor.
 * @param[in] TransA    flag indicating whether A should be transposed before multiplication
 * @param[in] TransB    flag indicating whether B should be transposed before multiplication
 * @param[in] M         number of rows in A/C
 * @param[in] N         number of columns in B/C
 * @param[in] K         number of rows in B, number of columns in A
 * @param[in] alpha     coefficient alpha
 * @param[in] A         matrix A
 * @param[in] lda       leading dimension of A (m or k)
 * @param[in] B         matrix B
 * @param[in] ldb       leading dimension of B (k or n)
 * @param[in] beta      coefficient beta
 * @param[in,out] C     matrix C
 * @param[in] ldc       leading dimension of C (m or n)
 * @return void
 */
void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    assert(layout == AMICI_BLAS_ColMajor);

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
 * amici_dgemm provides an interface to the CBlas matrix vector multiplication routine dgemv. This routines computes
 * y = alpha*A*x + beta*y with A: [MxN] x:[Nx1] y:[Mx1]
 *
 * @param[in] layout    always needs to be AMICI_BLAS_ColMajor.
 * @param[in] TransA    flag indicating whether A should be transposed before multiplication
 * @param[in] M         number of rows in A
 * @param[in] N         number of columns in A
 * @param[in] alpha     coefficient alpha
 * @param[in] A         matrix A
 * @param[in] lda       leading dimension of A (m or n)
 * @param[in] X         vector X
 * @param[in] incX      increment for entries of X
 * @param[in] beta      coefficient beta
 * @param[in,out] Y     vector Y
 * @param[in] incY      increment for entries of Y
 * @return void
 */
void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    assert(layout == AMICI_BLAS_ColMajor);

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
 * expDataFromMatlabCall parses the experimental data from the matlab call and writes it to an ExpData class object
 *
 * @param[in] prhs pointer to the array of input arguments @type mxArray
 * @param[in] udata pointer to user data object @type UserData
 * @param[in] model pointer to the model object, this is necessary to perform dimension checks @type Model
 * @return edata pointer to experimental data object @type ExpData
 */
ExpData *expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata, Model *model) {

    int nt_my = 0, ny_my = 0, nt_sigmay = 0,
        ny_sigmay = 0; /* integers with problem dimensionality */
    int ne_mz = 0, nz_mz = 0, ne_sigmaz = 0,
        nz_sigmaz = 0; /* integers with problem dimensionality */

    char errmsg[200];

    ExpData *edata = new ExpData(udata, model);
    if (edata->my == NULL) {
        goto freturn;
    }

    // Data provided / required?
    if (!mxGetPr(prhs[8])) {
        if (udata->sensi >= AMICI_SENSI_ORDER_FIRST &&
            udata->sensi_meth == AMICI_SENSI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data", "No data provided!");
            goto freturn;
        }
        goto freturn;
    }

    if (mxGetProperty(prhs[8], 0, "Y")) {
        ny_my = (int)mxGetN(mxGetProperty(prhs[8], 0, "Y"));
        if (ny_my != model->nytrue) {
            sprintf(errmsg, "Number of observables in data matrix (%i) does "
                            "not match model ny (%i)",
                    ny_my, model->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nyy", errmsg);
            goto freturn;
        }
        nt_my = (int)mxGetM(mxGetProperty(prhs[8], 0, "Y"));
        if (nt_my != udata->nt) {
            sprintf(errmsg, "Number of time-points in data matrix does (%i) "
                            "not match provided time vector (%i)",
                    nt_my, udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:nty", errmsg);
            goto freturn;
        }
        memcpy(edata->my, mxGetPr(mxGetProperty(prhs[8], 0, "Y")),
               ny_my * nt_my * sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Y",
                       "Field Y not specified as field in data struct!");
        goto freturn;
    }

    if (mxGetProperty(prhs[8], 0, "Sigma_Y")) {
        ny_sigmay = (int)mxGetN(mxGetProperty(prhs[8], 0, "Sigma_Y"));
        if (ny_sigmay != model->nytrue) {
            sprintf(errmsg, "Number of observables in data-sigma matrix (%i) "
                            "does not match model ny (%i)",
                    ny_sigmay, model->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nysdy", errmsg);
            goto freturn;
        }
        nt_sigmay = (int)mxGetM(mxGetProperty(prhs[8], 0, "Sigma_Y"));
        if (nt_sigmay != udata->nt) {
            sprintf(errmsg, "Number of time-points in data-sigma matrix (%i) "
                            "does not match provided time vector (%i)",
                    nt_sigmay, udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:ntsdy", errmsg);
            goto freturn;
        }
        memcpy(edata->sigmay, mxGetPr(mxGetProperty(prhs[8], 0, "Sigma_Y")),
               ny_sigmay * nt_sigmay * sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Y",
                       "Field Sigma_Y not specified as field in data struct!");
        goto freturn;
    }
    if (mxGetProperty(prhs[8], 0, "Z")) {
        nz_mz = (int)mxGetN(mxGetProperty(prhs[8], 0, "Z"));
        if (nz_mz != model->nztrue) {
            sprintf(errmsg, "Number of events in event matrix (%i) does not "
                            "match provided nz (%i)",
                    nz_mz, model->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nenz", errmsg);
            goto freturn;
        }
        ne_mz = (int)mxGetM(mxGetProperty(prhs[8], 0, "Z"));
        if (ne_mz != udata->nmaxevent) {
            sprintf(errmsg, "Number of time-points in event matrix (%i) does "
                            "not match provided nmaxevent (%i)",
                    ne_mz, udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnz", errmsg);
            goto freturn;
        }
        memcpy(edata->mz, mxGetPr(mxGetProperty(prhs[8], 0, "Z")),
               nz_mz * ne_mz * sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Z",
                       "Field Z not specified as field in data struct!");
        goto freturn;
    }

    if (mxGetProperty(prhs[8], 0, "Sigma_Z")) {
        nz_sigmaz = (int)mxGetN(mxGetProperty(prhs[8], 0, "Sigma_Z"));
        if (nz_sigmaz != model->nztrue) {
            sprintf(errmsg, "Number of events in event-sigma matrix (%i) does "
                            "not match provided nz (%i)",
                    nz_sigmaz, model->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nensdz", errmsg);
            goto freturn;
        }
        ne_sigmaz = (int)mxGetM(mxGetProperty(prhs[8], 0, "Sigma_Z"));
        if (ne_sigmaz != udata->nmaxevent) {
            sprintf(errmsg, "Number of time-points in event-sigma matrix (%i) "
                            "does not match provided nmaxevent (%i)",
                    ne_sigmaz, udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz", errmsg);
            goto freturn;
        }
        memcpy(edata->sigmaz, mxGetPr(mxGetProperty(prhs[8], 0, "Sigma_Z")),
               ne_sigmaz * nz_sigmaz * sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Z",
                       "Field Sigma_Z not specified as field in data struct!");
        goto freturn;
    }

    return edata;
freturn:
    if (edata)
        delete edata;
    return NULL;
}
