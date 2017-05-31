#include "include/amici_interface_matlab.h"

#include "wrapfunctions.h" /* user functions */

#include <cstring>
#include <assert.h>
#include <blas.h>

#include <include/edata_accessors.h>
#include <include/udata_accessors.h>
#include <include/rdata_accessors.h>
#include <include/tdata_accessors.h>

/**
 * @ brief initialise matrix and attach to the field
 * @ param FIELD name of the field to which the matrix will be attached
 * @ param D1 number of rows in the matrix
 * @ param D2 number of columns in the matrix
 */
#define initField2(FIELD,D1,D2) \
    mxArray *mx ## FIELD = mxCreateDoubleMatrix(D1,D2,mxREAL); \
    FIELD ## data = mxGetPr(mx ## FIELD); \
    mxSetField(mxsol,0,#FIELD,mx ## FIELD)

/**
 * @ brief initialise 3D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 */
#define initField3(FIELD,D1,D2,D3) \
    mwSize dims ## FIELD[] = {(mwSize)(D1), (mwSize)(D2), (mwSize)(D3)}; \
    mxArray *mx ## FIELD = mxCreateNumericArray(3,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
    FIELD ## data = mxGetPr(mx ## FIELD); \
    mxSetField(mxsol,0,#FIELD,mx ## FIELD)

/**
 * @ brief initialise 4D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 * @ param D4 number of elements in the fourth dimension of the tensor
 */
#define initField4(FIELD,D1,D2,D3,D4) \
    mwSize dims ## FIELD[] = {(mwSize)(D1), (mwSize)(D2), (mwSize)(D3), (mwSize)(D4)}; \
    mxArray *mx ## FIELD = mxCreateNumericArray(4,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
    FIELD ## data = mxGetPr(mx ## FIELD); \
    mxSetField(mxsol,0,#FIELD,mx ## FIELD)


/**
 * @ brief extract information from a property of a matlab class (scalar)
 * @ param OPTION name of the property
 * @ param TYPE class to which the information should be cast
 */
#define readOptionScalar(OPTION,TYPE) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3],0,#OPTION)); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
        return(NULL); \
    }


/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION) \
    if(mxGetProperty(prhs[3],0,#OPTION)){ \
        OPTION = (double *) mxGetData(mxGetProperty(prhs[3],0,#OPTION)); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
        return(NULL); \
    }



UserData *userDataFromMatlabCall(const mxArray *prhs[]) {
    /* User udata structure */
    UserData *udata = new UserData();
    if(udata==NULL) return NULL;

    init_modeldims(udata);

    /* time */

    if (!prhs[0]) {
        errMsgIdAndTxt("AMICI:mex:tout","No time vector provided!");
    }
    ts = mxGetPr(prhs[0]);

    nt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);

    /* parameters */

    if (!prhs[1]) {
        errMsgIdAndTxt("AMICI:mex:theta","No parameter vector provided!");
    }
    p = mxGetPr(prhs[1]);

    /* constants */

    if (!prhs[2]) {
        errMsgIdAndTxt("AMICI:mex:kappa","No constant vector provided!");
    }
    k = mxGetPr(prhs[2]);

    if (!prhs[3]) {
        errMsgIdAndTxt("AMICI:mex:options","No options provided!");
    }

    /* plist */
    realtype *plistdata;
    if (!prhs[4]) {
        errMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
    } else {
        nplist = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);
        plistdata = mxGetPr(prhs[4]);
    }

    plist = new int[nplist]();
    for (int ip=0; ip<nplist; ip++) {
        plist[ip] = (int)plistdata[ip];
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
        idlist = (double *) mxGetData(mxGetProperty(prhs[3],0,"id")); \
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
    }

    pbar = mxGetPr(prhs[5]);

    /* xscale */
    if (!prhs[6]) {
        errMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
    }

    xbar = mxGetPr(prhs[6]);

    /* Check, if initial states and sensitivities are passed by user or must be calculated */
    x0data = NULL;
    sx0data = NULL;
    if (prhs[7]) {
        if(mxGetField(prhs[7], 0 ,"x0")) {
            if ((mxGetM(mxGetField(prhs[7], 0 ,"x0")) * mxGetN(mxGetField(prhs[7], 0 ,"x0")))>0) {
                x0data = mxGetPr(mxGetField(prhs[7], 0 ,"x0"));

                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"x0")) != 1) { errMsgIdAndTxt("AMICI:mex:x0","Number of rows in x0 field must be equal to 1!"); }
                if(mxGetM(mxGetField(prhs[7], 0 ,"x0")) != nx) { errMsgIdAndTxt("AMICI:mex:x0","Number of columns in x0 field does not agree with number of model states!"); }
            }
        }

        if(mxGetField(prhs[7], 0 ,"sx0")) {
            if ((mxGetM(mxGetField(prhs[7], 0 ,"sx0")) * mxGetN(mxGetField(prhs[7], 0 ,"sx0")))>0) {
                sx0data = mxGetPr(mxGetField(prhs[7], 0 ,"sx0"));

                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"sx0")) != nplist) { errMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); }
                if(mxGetM(mxGetField(prhs[7], 0 ,"sx0")) != nx) { errMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); }
            }
        }
    }

    processUserData(udata);

    return(udata);
}


ReturnData *setupReturnData(mxArray *plhs[], const UserData *udata, double *pstatus) {

    const char *field_names_sol[] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};

    /* Return rdata structure */
    ReturnData *rdata = (ReturnData*) mxMalloc(sizeof *rdata);
    if (rdata == NULL) return(NULL);

    memset(rdata, 0, sizeof(*rdata));

    mxArray *mxsol = mxCreateStructMatrix(1,1,29,field_names_sol);
    plhs[0] = mxsol;

    #include "include/amici_init_return_data_fields.h"

    mxArray *mxstatus = mxCreateDoubleMatrix(1,1,mxREAL);
    mxSetPr(mxstatus,pstatus);
    mxSetField(mxsol,0,"status",mxstatus);

    mxArray *mxts = mxCreateDoubleMatrix(nt,1,mxREAL);
    tsdata = mxGetPr(mxts);
    mxSetField(mxsol,0,"t",mxts);

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
        if(sensi >= AMI_SENSI_ORDER_FIRST && sensi_meth == AMI_SENSI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data","No data provided!");
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

        if (nmyt != nt) {
            sprintf(errmsg,"Number of time-points in data matrix does (%i) not match provided time vector (%i)",nmyt,nt);
            errMsgIdAndTxt("AMICI:mex:data:nty",errmsg);
            return NULL;
        }

        if (nysigmat != nt) {
            sprintf(errmsg,"Number of time-points in data-sigma matrix (%i) does not match provided time vector (%i)",nysigmat,nt);
            errMsgIdAndTxt("AMICI:mex:data:ntsdy",errmsg);
            return NULL;
        }

        if (nmyy != nytrue) {
            sprintf(errmsg,"Number of observables in data matrix (%i) does not match model ny (%i)",nmyy,nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nyy",errmsg);
            return NULL;
        }

        if (nysigmay != nytrue) {
            sprintf(errmsg,"Number of observables in data-sigma matrix (%i) does not match model ny (%i)",nysigmay,nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nysdy",errmsg);
            return NULL;
        }

        if (nmzt != nmaxevent) {
            sprintf(errmsg,"Number of time-points in event matrix (%i) does not match provided nmaxevent (%i)",nmzt,nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnz",errmsg);
            return NULL;
        }

        if (nzsigmat != nmaxevent) {
            sprintf(errmsg,"Number of time-points in event-sigma matrix (%i) does not match provided nmaxevent (%i)",nzsigmat,nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz",errmsg);
            return NULL;
        }

        if (nmzy != nztrue) {
            sprintf(errmsg,"Number of events in event matrix (%i) does not match provided nz (%i)",nmzy,nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nenz",errmsg);
            return NULL;
        }

        if (nzsigmay != nztrue) {
            sprintf(errmsg,"Number of events in event-sigma matrix (%i) does not match provided nz (%i)",nzsigmay,nztrue);
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

    const long int M_ = M;
    const long int N_ = N;
    const long int K_ = K;
    const long int lda_ = lda;
    const long int ldb_ = ldb;
    const long int ldc_ = ldc;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);
    const char transB = amici_blasCBlasTransToBlasTrans(TransB);

#if defined(_WIN32)
    dgemm(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_);
#else
    dgemm_(&transA, &transB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_);
#endif
}

void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta, double *Y, const int incY)
{
    assert(layout == AMICI_BLAS_ColMajor);

    const long int M_ = M;
    const long int N_ = N;
    const long int lda_ = lda;
    const long int incX_ = incX;
    const long int incY_ = incY;
    const char transA = amici_blasCBlasTransToBlasTrans(TransA);

    assert(layout == AMICI_BLAS_ColMajor);
#if defined(_WIN32)
    dgemv(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#else
    dgemv_(&transA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_);
#endif
}
