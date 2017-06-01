#include "include/amici_interface_matlab.h"

#include "wrapfunctions.h" /* user functions */

#include <cstring>
#include <assert.h>
#include <blas.h>
#include <new>
#include <include/edata_accessors.h>
#include <include/tdata_accessors.h>

/**
 * @ brief extract information from a property of a matlab class (scalar)
 * @ param OPTION name of the property
 * @ param TYPE class to which the information should be cast
 */
#define readOptionScalar(OPTION,TYPE) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        udata->OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3],0,#OPTION)); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
        return(NULL); \
    }


/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        udata->OPTION = (double *) mxGetData(mxGetProperty(prhs[3],0,#OPTION)); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
        return(NULL); \
    }



UserData *userDataFromMatlabCall(const mxArray *prhs[]) {
    /* User udata structure */
    UserData *udata = getUserData();
    if(udata==NULL) return NULL;

    /* time */

    if (!prhs[0]) {
        errMsgIdAndTxt("AMICI:mex:tout","No time vector provided!");
        return NULL;
    }
    udata->ts = mxGetPr(prhs[0]);

    udata->nt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);

    /* parameters */

    if (!prhs[1]) {
        errMsgIdAndTxt("AMICI:mex:theta","No parameter vector provided!");
        return NULL;
    }
    udata->p = mxGetPr(prhs[1]);

    /* constants */

    if (!prhs[2]) {
        errMsgIdAndTxt("AMICI:mex:kappa","No constant vector provided!");
        return NULL;
    }
    udata->k = mxGetPr(prhs[2]);

    if (!prhs[3]) {
        errMsgIdAndTxt("AMICI:mex:options","No options provided!");
        return NULL;
    }

    /* plist */
    realtype *plistdata;
    if (!prhs[4]) {
        errMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
        return NULL;
    } else {
        udata->nplist = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);
        plistdata = mxGetPr(prhs[4]);
    }

    udata->plist = new int[udata->nplist]();
    for (int ip=0; ip<udata->nplist; ip++) {
        udata->plist[ip] = (int)plistdata[ip];
    }

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

    if(mxGetProperty(prhs[3],0,"id")){ \
        udata->idlist = (double *) mxGetData(mxGetProperty(prhs[3],0,"id")); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
        return(NULL); \
    }

    readOptionData(z2event)
    readOptionData(qpositivex)
    readOptionScalar(sensi,AMI_sensi_order)
    readOptionScalar(ism,int)
    readOptionScalar(sensi_meth,AMI_sensi_meth)
    readOptionScalar(ordering,int)

    /* pbar */
    if (!prhs[5]) {
        errMsgIdAndTxt("AMICI:mex:pbar","No parameter scales provided!");
        return NULL;
    }

    udata->pbar = mxGetPr(prhs[5]);

    /* xscale */
    if (!prhs[6]) {
        errMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
        return NULL;
    }

    udata->xbar = mxGetPr(prhs[6]);

    /* Check, if initial states and sensitivities are passed by user or must be calculated */
    if (prhs[7]) {
        if(mxGetField(prhs[7], 0 ,"x0")) {
            if ((mxGetM(mxGetField(prhs[7], 0 ,"x0")) * mxGetN(mxGetField(prhs[7], 0 ,"x0")))>0) {
                udata->x0data = mxGetPr(mxGetField(prhs[7], 0 ,"x0"));

                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"x0")) != 1) { errMsgIdAndTxt("AMICI:mex:x0","Number of rows in x0 field must be equal to 1!"); return NULL; }
                if(mxGetM(mxGetField(prhs[7], 0 ,"x0")) != udata->nx) { errMsgIdAndTxt("AMICI:mex:x0","Number of columns in x0 field does not agree with number of model states!"); return NULL; }
            }
        }

        if(mxGetField(prhs[7], 0 ,"sx0")) {
            if ((mxGetM(mxGetField(prhs[7], 0 ,"sx0")) * mxGetN(mxGetField(prhs[7], 0 ,"sx0")))>0) {
                udata->sx0data = mxGetPr(mxGetField(prhs[7], 0 ,"sx0"));

                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"sx0")) != udata->nplist) { errMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); return NULL; }
                if(mxGetM(mxGetField(prhs[7], 0 ,"sx0")) != udata->nx) { errMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); return NULL; }
            }
        }
    }

    return(udata);
}


ReturnDataMatlab *setupReturnData(mxArray *plhs[], const UserData *udata, double *pstatus) {

    ReturnDataMatlab *rdata = (ReturnDataMatlab*) mxMalloc(sizeof *rdata);
    if (rdata == NULL) return(NULL);
    new(rdata) ReturnDataMatlab(udata);
    plhs[0] = rdata->mxsol;

    return(rdata);
}


ExpData *expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata, int *status) {
    int nmyt = 0, nmyy = 0, nysigmat = 0, nysigmay = 0; /* integers with problem dimensionality */
    int nmzt = 0, nmzy = 0, nzsigmat = 0, nzsigmay = 0; /* integers with problem dimensionality */

    char *errmsg;
    errmsg = new char[200]();

    *status = -97;

    ExpData *edata; /* returned rdata struct */

    if ((mxGetM(prhs[8]) == 0 && mxGetN(prhs[8]) == 0) || !prhs[8]) {
        if(udata->sensi >= AMI_SENSI_ORDER_FIRST && udata->sensi_meth == AMI_SENSI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data","No data provided!");
            return NULL;
        } else {
            *status = 0;
        }
        return NULL;
    } else {
        /* Return edata structure */
        edata = new ExpData();
        if (edata == NULL) return(NULL);

        if (mxGetProperty(prhs[8], 0 ,"Y")) {
            my = mxGetPr(mxGetProperty(prhs[8], 0 ,"Y"));
            nmyy = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Y"));
            nmyt = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Y"));
        } else {
            errMsgIdAndTxt("AMICI:mex:data:Y","Field Y not specified as field in data struct!");
            return NULL;
        }

        if (mxGetProperty(prhs[8], 0 ,"Sigma_Y")) {
            ysigma = mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
            nysigmay = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
            nysigmat = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
        } else {
            errMsgIdAndTxt("AMICI:mex:data:Sigma_Y","Field Sigma_Y not specified as field in data struct!");
            return NULL;
        }
        if (mxGetProperty(prhs[8], 0 ,"Z")) {
            mz = mxGetPr(mxGetProperty(prhs[8], 0 ,"Z"));
            nmzy = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Z"));
            nmzt = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Z"));
        } else {
            errMsgIdAndTxt("AMICI:mex:data:Z","Field Z not specified as field in data struct!");
            return NULL;
        }

        if (mxGetProperty(prhs[8], 0 ,"Sigma_Z")) {
            zsigma = mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
            nzsigmay = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
            nzsigmat = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
        } else {
            errMsgIdAndTxt("AMICI:mex:data:Sigma_Z","Field Sigma_Z not specified as field in data struct!");
            return NULL;
        }

        if (nmyt != udata->nt) {
            sprintf(errmsg,"Number of time-points in data matrix does (%i) not match provided time vector (%i)",nmyt,udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:nty",errmsg);
            return NULL;
        }

        if (nysigmat != udata->nt) {
            sprintf(errmsg,"Number of time-points in data-sigma matrix (%i) does not match provided time vector (%i)",nysigmat,udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:ntsdy",errmsg);
            return NULL;
        }

        if (nmyy != udata->nytrue) {
            sprintf(errmsg,"Number of observables in data matrix (%i) does not match model ny (%i)",nmyy,udata->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nyy",errmsg);
            return NULL;
        }

        if (nysigmay != udata->nytrue) {
            sprintf(errmsg,"Number of observables in data-sigma matrix (%i) does not match model ny (%i)",nysigmay,udata->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nysdy",errmsg);
            return NULL;
        }

        if (nmzt != udata->nmaxevent) {
            sprintf(errmsg,"Number of time-points in event matrix (%i) does not match provided nmaxevent (%i)",nmzt,udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnz",errmsg);
            return NULL;
        }

        if (nzsigmat != udata->nmaxevent) {
            sprintf(errmsg,"Number of time-points in event-sigma matrix (%i) does not match provided nmaxevent (%i)",nzsigmat,udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz",errmsg);
            return NULL;
        }

        if (nmzy != udata->nztrue) {
            sprintf(errmsg,"Number of events in event matrix (%i) does not match provided nz (%i)",nmzy,udata->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nenz",errmsg);
            return NULL;
        }

        if (nzsigmay != udata->nztrue) {
            sprintf(errmsg,"Number of events in event-sigma matrix (%i) does not match provided nz (%i)",nzsigmay,udata->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nensdz",errmsg);
            return NULL;
        }
    }

    delete[] errmsg;

    *status = 0;

    return(edata);
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

ReturnDataMatlab::ReturnDataMatlab(const UserData *udata) : ReturnData(udata)
{

}

void ReturnDataMatlab::initFields(const UserData *udata)
{
    const char *field_names_sol[] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};
    mxsol = mxCreateStructMatrix(1,1,29,field_names_sol);

    ReturnData::initFields(udata);
}

void ReturnDataMatlab::initField1(double **fieldPointer, const char *fieldName, int dim)
{
    mxArray *array = mxCreateDoubleMatrix(dim, 1, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

}

void ReturnDataMatlab::initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2)
{
    mxArray *array = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);
}

void ReturnDataMatlab::initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    mxArray *array = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

}

void ReturnDataMatlab::initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3), (mwSize)(dim4)};
    mxArray *array = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
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
