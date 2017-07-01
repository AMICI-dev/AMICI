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
udata->OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3],0,#OPTION)); \
} else { \
warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
return NULL; \
}

/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION) \
if(mxGetProperty(prhs[3],0,#OPTION)){ \
mxArray *a = mxGetProperty(prhs[3],0,#OPTION); \
int len = (int) mxGetM(a) * mxGetN(a); \
udata->OPTION = new double[len]; \
memcpy(udata->OPTION, mxGetData(a), sizeof(double) * len); \
} else { \
warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options do not have field " #OPTION "!"); \
return NULL; \
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

    UserData *udata = userDataFromMatlabCall(prhs);

    ReturnDataMatlab *rdata = new ReturnDataMatlab(udata);
    plhs[0] = rdata->mxsol;
    if(*(rdata->status) != AMICI_SUCCESS)
        return;

    ExpData *edata = expDataFromMatlabCall(prhs, udata);
    
    if(udata->nx > 0) {
        *(rdata->status) = (double) runAmiciSimulation(udata, edata, rdata);
    }

    if(udata) delete udata;
    if(rdata) delete rdata;
    if(edata) delete edata;
}

UserData *userDataFromMatlabCall(const mxArray *prhs[]) {
    
    UserData *udata = new UserData(getUserData());
    
    /* time */
    if (prhs[0]) {
        udata->nt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);
        udata->ts = new realtype[udata->nt];
        memcpy(udata->ts, mxGetPr(prhs[0]), sizeof(realtype) * udata->nt);
    } else {
        errMsgIdAndTxt("AMICI:mex:tout","No time vector provided!");
        return NULL;
    }
    
    /* parameters */
    if (prhs[1]) {
        udata->p = new realtype[udata->np];
        memcpy(udata->p, mxGetPr(prhs[1]), sizeof(realtype) * udata->np);
    } else {
        errMsgIdAndTxt("AMICI:mex:theta","No parameter vector provided!");
        return NULL;
    }
    
    /* constants */
    if (prhs[2]) {
        int lenK = (int) mxGetM(prhs[2]) * mxGetN(prhs[2]);
        udata->k = new realtype[lenK];
        memcpy(udata->k, mxGetPr(prhs[2]), sizeof(realtype) * lenK);
        
    } else {
        errMsgIdAndTxt("AMICI:mex:kappa","No constant vector provided!");
        return NULL;
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
            udata->idlist = new realtype[lenIdlist];
            memcpy(udata->idlist, mxGetData(idlist), sizeof(realtype) * lenIdlist);
        } else {
            warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!");
            return NULL;
        }
        
        readOptionData(z2event)
        readOptionData(qpositivex)
        readOptionScalar(sensi,AMICI_sensi_order)
        readOptionScalar(ism,int)
        readOptionScalar(sensi_meth,AMICI_sensi_meth)
        readOptionScalar(ordering,int)
        
    } else {
        errMsgIdAndTxt("AMICI:mex:options","No options provided!");
        return NULL;
    }
    
    /* plist */
    if (prhs[4]) {
        udata->nplist = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);
        udata->plist = new int[udata->nplist]();
        realtype *plistdata = mxGetPr(prhs[4]);
        
        for (int ip = 0; ip < udata->nplist; ip++) {
            udata->plist[ip] = (int)plistdata[ip];
        }
        
    } else {
        errMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
        return NULL;
    }
    
    /* pbar */
    if (prhs[5]) {
        int lenPBar = (int) mxGetM(prhs[5]) * mxGetN(prhs[5]);
        udata->pbar = new realtype[lenPBar]();
        memcpy(udata->pbar, mxGetPr(prhs[5]), sizeof(realtype) * lenPBar);
    } else {
        errMsgIdAndTxt("AMICI:mex:pbar","No parameter scales provided!");
        return NULL;
    }
    
    /* xscale */
    if (prhs[6]) {
        int lenXBar = (int) mxGetM(prhs[6]) * mxGetN(prhs[6]);
        udata->xbar = new realtype[lenXBar]();
        memcpy(udata->xbar, mxGetPr(prhs[6]), sizeof(realtype) * lenXBar);
    } else {
        errMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
        return NULL;
    }
    
    /* Check, if initial states and sensitivities are passed by user or must be calculated */
    if (prhs[7]) {
        mxArray *x0 = mxGetField(prhs[7], 0 ,"x0");
        if(x0) {
            /* check dimensions */
            if(mxGetN(x0) != 1) { errMsgIdAndTxt("AMICI:mex:x0","Number of rows in x0 field must be equal to 1!"); return NULL; }
            if(mxGetM(x0) != udata->nx) { errMsgIdAndTxt("AMICI:mex:x0","Number of columns in x0 field does not agree with number of model states!"); return NULL; }
            
            if ((mxGetM(x0) * mxGetN(x0)) > 0) {
                udata->x0data = new realtype[udata->nx];
                memcpy(udata->x0data, mxGetPr(x0), sizeof(realtype) * udata->nx);
            }
        }
        
        mxArray *sx0 = mxGetField(prhs[7], 0 ,"sx0");
        if(sx0 && (mxGetM(sx0) * mxGetN(sx0)) > 0) {
            /* check dimensions */
            if(mxGetN(sx0) != udata->nplist) { errMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); return NULL; }
            if(mxGetM(sx0) != udata->nx) { errMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); return NULL; }
            
            udata->sx0data = new realtype[udata->nx * udata->nplist];
            memcpy(udata->sx0data, mxGetPr(sx0), sizeof(realtype) * udata->nx * udata->nplist);
        }
    }
    return udata;
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
    const int numFields = 29;
    const char *field_names_sol[numFields] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};

    mxsol = mxCreateStructMatrix(1, 1, numFields, field_names_sol);

    ReturnData::initFields(udata);
}

void ReturnDataMatlab::initField1(double **fieldPointer, const char *fieldName, int dim)
{
    mxArray *array = mxCreateDoubleMatrix(dim, 1, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if(status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2)
{
    mxArray *array = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if(status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    mxArray *array = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if(status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4)
{
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3), (mwSize)(dim4)};
    mxArray *array = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if(status && array == NULL)
        *status = AMICI_ERROR_RDATA;
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

ExpData *expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata) {
    
    ExpData *edata = new ExpData(udata);
    
    // Data provided / required?
    if ((!prhs[8] || (mxGetM(prhs[8]) == 0 && mxGetN(prhs[8]) == 0))) {
        if(udata->sensi >= AMICI_SENSI_ORDER_FIRST && udata->sensi_meth == AMICI_SENSI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data","No data provided!");
            return NULL;
        }
        return NULL;
    }
    
    int nt_my = 0, ny_my = 0, nt_sigmay = 0, ny_sigmay = 0; /* integers with problem dimensionality */
    int ne_mz = 0, nz_mz = 0, ne_sigmaz = 0, nz_sigmaz = 0; /* integers with problem dimensionality */
    
    char errmsg[200];
    
    if (mxGetProperty(prhs[8], 0 ,"Y")) {
        ny_my = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Y"));
        if (ny_my != udata->nytrue) {
            sprintf(errmsg,"Number of observables in data matrix (%i) does not match model ny (%i)",ny_my,udata->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nyy", errmsg);
            return NULL;
        }
        nt_my = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Y"));
        if (nt_my != udata->nt) {
            sprintf(errmsg,"Number of time-points in data matrix does (%i) not match provided time vector (%i)",nt_my,udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:nty", errmsg);
            return NULL;
        }
        memcpy(edata->my,mxGetPr(mxGetProperty(prhs[8], 0 ,"Y")),ny_my*nt_my*sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Y","Field Y not specified as field in data struct!");
        return NULL;
    }
    
    if (mxGetProperty(prhs[8], 0 ,"Sigma_Y")) {
        ny_sigmay = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
        if (ny_sigmay != udata->nytrue) {
            sprintf(errmsg,"Number of observables in data-sigma matrix (%i) does not match model ny (%i)",ny_sigmay,udata->nytrue);
            errMsgIdAndTxt("AMICI:mex:data:nysdy", errmsg);
            return NULL;
        }
        nt_sigmay = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Y"));
        if (nt_sigmay != udata->nt) {
            sprintf(errmsg,"Number of time-points in data-sigma matrix (%i) does not match provided time vector (%i)",nt_sigmay,udata->nt);
            errMsgIdAndTxt("AMICI:mex:data:ntsdy", errmsg);
            return NULL;
        }
        memcpy(edata->sigmay,mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Y")),ny_sigmay*nt_sigmay*sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Y","Field Sigma_Y not specified as field in data struct!");
        return NULL;
    }
    if (mxGetProperty(prhs[8], 0 ,"Z")) {
        nz_mz = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Z"));
        if (nz_mz != udata->nztrue) {
            sprintf(errmsg,"Number of events in event matrix (%i) does not match provided nz (%i)",nz_mz,udata->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nenz", errmsg);
            return NULL;
        }
        ne_mz = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Z"));
        if (ne_mz != udata->nmaxevent) {
            sprintf(errmsg,"Number of time-points in event matrix (%i) does not match provided nmaxevent (%i)",ne_mz,udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnz", errmsg);
            return NULL;
        }
        memcpy(edata->mz,mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Y")),nz_mz*ne_mz*sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Z","Field Z not specified as field in data struct!");
        return NULL;
    }
    
    if (mxGetProperty(prhs[8], 0 ,"Sigma_Z")) {
        nz_sigmaz = (int) mxGetN(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
        if (nz_sigmaz != udata->nztrue) {
            sprintf(errmsg,"Number of events in event-sigma matrix (%i) does not match provided nz (%i)",nz_sigmaz,udata->nztrue);
            errMsgIdAndTxt("AMICI:mex:data:nensdz", errmsg);
            return NULL;
        }
        ne_sigmaz = (int) mxGetM(mxGetProperty(prhs[8], 0 ,"Sigma_Z"));
        if (ne_sigmaz != udata->nmaxevent) {
            sprintf(errmsg,"Number of time-points in event-sigma matrix (%i) does not match provided nmaxevent (%i)",ne_sigmaz,udata->nmaxevent);
            errMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz", errmsg);
            return NULL;
        }
        memcpy(edata->sigmaz,mxGetPr(mxGetProperty(prhs[8], 0 ,"Sigma_Z")),ne_sigmaz*nz_sigmaz*sizeof(double));
    } else {
        errMsgIdAndTxt("AMICI:mex:data:Sigma_Z","Field Sigma_Z not specified as field in data struct!");
        return NULL;
    }
    return edata;
}

