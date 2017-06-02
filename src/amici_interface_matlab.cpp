/**
 * @file   amiwrap.cpp
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the mex file from matlab
 */

#include "include/amici_interface_matlab.h"
#include "wrapfunctions.h" /* user functions */

#include <cstring>
#include <assert.h>
#include <blas.h>
#include <include/edata_accessors.h>
#include <include/tdata_accessors.h>

#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

/**
 * @ brief extract information from a property of a matlab class (scalar)
 * @ param OPTION name of the property
 * @ param TYPE class to which the information should be cast
 */
#define readOptionScalar(OPTION,TYPE) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        udata.OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3],0,#OPTION)); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
        return(udata); \
    }


/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        mxArray *a = mxGetProperty(prhs[3],0,#OPTION); \
        int len = (int) mxGetM(a) * mxGetN(a); \
        udata.OPTION = new realtype[len]; \
        memcpy(udata.OPTION, mxGetData(a), sizeof(realtype) * len); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
        return(udata); \
    }



/*!
 * mexFunction is the main function of the mex simulation file this function carries out all numerical integration and writes results into the sol struct.
 *
 * @param[in] nlhs number of output arguments of the matlab call @type int
 * @param[out] plhs pointer to the array of output arguments @type mxArray
 * @param[in] nrhs number of input arguments of the matlab call @type int
 * @param[in] prhs pointer to the array of input arguments @type mxArray
 * @return void
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* return status flag */
    int status = 0;

    UserData udata = userDataFromMatlabCall(prhs, &status);
    ReturnDataMatlab rdata(&udata);

    plhs[0] = rdata.mxsol;

    if(status == 0 && udata.nx > 0) {
        ExpData edata = expDataFromMatlabCall(prhs, &udata, &status);
        if (status == 0)
            runAmiciSimulation(&udata, &edata, &rdata, &status);
    }

    *rdata.status = (double) status;
}


UserData userDataFromMatlabCall(const mxArray *prhs[], int *status) {
    *status = -98;

    UserData udata = getUserData();

    /* time */
    if (prhs[0]) {
        udata.nt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);
        udata.ts = new realtype[udata.nt];
        memcpy(udata.ts, mxGetPr(prhs[0]), sizeof(realtype) * udata.nt);
    } else {
        errMsgIdAndTxt("AMICI:mex:tout","No time vector provided!");
        return udata;
    }

    /* parameters */
    if (prhs[1]) {
        udata.p = new realtype[udata.np];
        memcpy(udata.p, mxGetPr(prhs[1]), sizeof(realtype) * udata.np);
    } else {
        errMsgIdAndTxt("AMICI:mex:theta","No parameter vector provided!");
        return udata;
    }

    /* constants */
    if (prhs[2]) {
        int lenK = (int) mxGetM(prhs[2]) * mxGetN(prhs[2]);
        udata.k = new realtype[lenK];
        memcpy(udata.k, mxGetPr(prhs[2]), sizeof(realtype) * lenK);

    } else {
        errMsgIdAndTxt("AMICI:mex:kappa","No constant vector provided!");
        return udata;
    }

    /* options */
    if (prhs[3]) {
        readOptionScalar(nmaxevent,int)
        readOptionScalar(tstart,double)
        readOptionScalar(atol,double)
        readOptionScalar(rtol,double)
        readOptionScalar(maxsteps,int)
        readOptionScalar(lmm,int)
        readOptionScalar(iter,int)
        readOptionScalar(interpType,int)
        readOptionScalar(linsol,int)
        readOptionScalar(stldet,booleantype)

        if(mxGetProperty(prhs[3],0,"id")){
            mxArray *idlist = mxGetProperty(prhs[3],0,"id");
            int lenIdlist = (int) mxGetM(idlist) * mxGetN(idlist);
            udata.idlist = new realtype[lenIdlist];
            memcpy(udata.idlist, mxGetData(idlist), sizeof(realtype) * lenIdlist);
        } else {
            warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!");
            return udata;
        }

        readOptionData(z2event)
        readOptionData(qpositivex)
        readOptionScalar(sensi,AMI_sensi_order)
        readOptionScalar(ism,int)
        readOptionScalar(sensi_meth,AMI_sensi_meth)
        readOptionScalar(ordering,int)

    } else {
        errMsgIdAndTxt("AMICI:mex:options","No options provided!");
        return udata;
    }

    /* plist */
    if (prhs[4]) {
        udata.nplist = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);
        udata.plist = new int[udata.nplist]();
        realtype *plistdata = mxGetPr(prhs[4]);

        for (int ip = 0; ip < udata.nplist; ip++) {
            udata.plist[ip] = (int)plistdata[ip];
        }

    } else {
        errMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
        return udata;
    }

    /* pbar */
    if (prhs[5]) {
        int lenPBar = (int) mxGetM(prhs[5]) * mxGetN(prhs[5]);
        udata.pbar = new realtype[lenPBar]();
        memcpy(udata.pbar, mxGetPr(prhs[5]), sizeof(realtype) * lenPBar);
    } else {
        errMsgIdAndTxt("AMICI:mex:pbar","No parameter scales provided!");
        return udata;
    }

    /* xscale */
    if (prhs[6]) {
        int lenXBar = (int) mxGetM(prhs[6]) * mxGetN(prhs[6]);
        udata.xbar = new realtype[lenXBar]();
        memcpy(udata.xbar, mxGetPr(prhs[6]), sizeof(realtype) * lenXBar);
    } else {
        errMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
        return udata;
    }

    /* Check, if initial states and sensitivities are passed by user or must be calculated */
    if (prhs[7]) {
        mxArray *x0 = mxGetField(prhs[7], 0 ,"x0");
        if(x0) {
            /* check dimensions */
            if(mxGetN(x0) != 1) { errMsgIdAndTxt("AMICI:mex:x0","Number of rows in x0 field must be equal to 1!"); return udata; }
            if(mxGetM(x0) != udata.nx) { errMsgIdAndTxt("AMICI:mex:x0","Number of columns in x0 field does not agree with number of model states!"); return udata; }

            if ((mxGetM(x0) * mxGetN(x0)) > 0) {
                udata.x0data = new realtype[udata.nx];
                memcpy(udata.x0data, mxGetPr(x0), sizeof(realtype) * udata.nx);
            }
        }

        mxArray *sx0 = mxGetField(prhs[7], 0 ,"sx0");
        if(sx0 && (mxGetM(sx0) * mxGetN(sx0)) > 0) {
            /* check dimensions */
            if(mxGetN(sx0) != udata.nplist) { errMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); return udata; }
            if(mxGetM(sx0) != udata.nx) { errMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); return udata; }

            udata.sx0data = new realtype[udata.nx * udata.nplist];
            memcpy(udata.sx0data, mxGetPr(sx0), sizeof(realtype) * udata.nx * udata.nplist);
        }
    }

    *status = 0;

    return(udata);
}

ExpData expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata, int *status) {
    *status = -97;
    ExpData _edata = {0};
    ExpData *edata = &_edata; // for accessor macros

    // Data provided / required?
    if ((!prhs[8] || (mxGetM(prhs[8]) == 0 && mxGetN(prhs[8]) == 0))) {
        if(udata->sensi >= AMI_SENSI_ORDER_FIRST && udata->sensi_meth == AMI_SENSI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data","No data provided!");
        } else {
            *status = 0;
        }
        return _edata;
    }

    int nmyt = 0, nmyy = 0, nysigmat = 0, nysigmay = 0; /* integers with problem dimensionality */
    int nmzt = 0, nmzy = 0, nzsigmat = 0, nzsigmay = 0; /* integers with problem dimensionality */

    if (mxGetProperty(prhs[8], 0 ,"Y")) {
        my = mxGetPr(mxGetProperty(prhs[8], 0 ,"Y"));
        nmyy = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Y"));
        nmyt = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Y"));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Y","Field Y not specified as field in data struct!");
        return _edata;
    }

    if (mxGetProperty(prhs[8], 0 ,"Sigma_Y")) {
        ysigma = mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
        nysigmay = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
        nysigmat = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Y","Field Sigma_Y not specified as field in data struct!");
        return _edata;
    }
    if (mxGetProperty(prhs[8], 0 ,"Z")) {
        mz = mxGetPr(mxGetProperty(prhs[8], 0 ,"Z"));
        nmzy = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Z"));
        nmzt = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Z"));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Z","Field Z not specified as field in data struct!");
        return _edata;
    }

    if (mxGetProperty(prhs[8], 0 ,"Sigma_Z")) {
        zsigma = mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
        nzsigmay = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
        nzsigmat = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Z","Field Sigma_Z not specified as field in data struct!");
        return _edata;
    }

    char errmsg[200];

    if (nmyt != udata->nt) {
        sprintf(errmsg,"Number of time-points in data matrix does (%i) not match provided time vector (%i)",nmyt,udata->nt);
        errMsgIdAndTxt("AMICI:mex:data:nty", errmsg);
    }

    if (nysigmat != udata->nt) {
        sprintf(errmsg,"Number of time-points in data-sigma matrix (%i) does not match provided time vector (%i)",nysigmat,udata->nt);
        errMsgIdAndTxt("AMICI:mex:data:ntsdy", errmsg);
    }

    if (nmyy != udata->nytrue) {
        sprintf(errmsg,"Number of observables in data matrix (%i) does not match model ny (%i)",nmyy,udata->nytrue);
        errMsgIdAndTxt("AMICI:mex:data:nyy", errmsg);
    }

    if (nysigmay != udata->nytrue) {
        sprintf(errmsg,"Number of observables in data-sigma matrix (%i) does not match model ny (%i)",nysigmay,udata->nytrue);
        errMsgIdAndTxt("AMICI:mex:data:nysdy", errmsg);
    }

    if (nmzt != udata->nmaxevent) {
        sprintf(errmsg,"Number of time-points in event matrix (%i) does not match provided nmaxevent (%i)",nmzt,udata->nmaxevent);
        errMsgIdAndTxt("AMICI:mex:data:nmaxeventnz", errmsg);
    }

    if (nzsigmat != udata->nmaxevent) {
        sprintf(errmsg,"Number of time-points in event-sigma matrix (%i) does not match provided nmaxevent (%i)",nzsigmat,udata->nmaxevent);
        errMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz", errmsg);
    }

    if (nmzy != udata->nztrue) {
        sprintf(errmsg,"Number of events in event matrix (%i) does not match provided nz (%i)",nmzy,udata->nztrue);
        errMsgIdAndTxt("AMICI:mex:data:nenz", errmsg);
    }

    if (nzsigmay != udata->nztrue) {
        sprintf(errmsg,"Number of events in event-sigma matrix (%i) does not match provided nz (%i)",nzsigmay,udata->nztrue);
        errMsgIdAndTxt("AMICI:mex:data:nensdz", errmsg);
    }

    *status = 0;

    return _edata;
}

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

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA, AMICI_BLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
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
    dgemm(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_);
#else
    dgemm_(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_);
#endif
}

ReturnDataMatlab::ReturnDataMatlab(const UserData *udata) : ReturnData()
{
    mxsol = NULL;
    freeFieldsOnDestruction = false;
    initFields(udata);
}

void ReturnDataMatlab::initFields(const UserData *udata)
{
    const char *field_names_sol[] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};
    mxsol = mxCreateStructMatrix(1,1,29,field_names_sol);

    ReturnData::initFields(udata);
}

void ReturnDataMatlab::initField1(double **fieldPointer, const char *fieldName, int dim)
{
    mxArray *array;
    array = mxCreateDoubleMatrix(dim, 1, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

}

void ReturnDataMatlab::initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2)
{
    mxArray *array;
    array = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);
}

void ReturnDataMatlab::initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    mxArray *array;
    array = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

}

void ReturnDataMatlab::initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3), (mwSize)(dim4)};
    mxArray *array;
    array = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

}

void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta, double *Y, const int incY)
{
    assert(layout == AMICI_BLAS_ColMajor);

    const ptrdiff_t M_ = M;
    const ptrdiff_t N_ = N;
    const ptrdiff_t lda_ = lda;
    const ptrdiff_t incX_ = incX;
    const ptrdiff_t incY_ = incY;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);

    assert(layout == AMICI_BLAS_ColMajor);
#if defined(_WIN32)
    dgemv(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#else
    dgemv_(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#endif
}
