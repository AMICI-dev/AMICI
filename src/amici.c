/**
 * @file   amici.c
 * @brief  core routines for integration
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <math.h>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif
#include <mex.h>
#include "wrapfunctions.h" /* user functions */
#include <include/amici.h> /* amici functions */

/** 
 * @ brief initialise matrix and attach to the field 
 * @ param FIELD name of the field to which the matrix will be attached
 * @ param D1 number of rows in the matrix
 * @ param D2 number of columns in the matrix
 */
#define initField2(FIELD,D1,D2) \
mx ## FIELD = mxCreateDoubleMatrix(D1,D2,mxREAL); \
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
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
mx ## FIELD = mxCreateNumericArray(3,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
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
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
dims ## FIELD[3]=D4; \
mx ## FIELD = mxCreateNumericArray(4,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
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
    mexWarnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
    return(NULL); \
}

/** 
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#define readOptionData(OPTION) \
if(mxGetProperty(prhs[3],0,#OPTION)){ \
    OPTION = mxGetData(mxGetProperty(prhs[3],0,#OPTION)); \
} else { \
    mexWarnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
    return(NULL); \
}

/** return value for successful execution */
#define AMI_SUCCESS               0

UserData setupUserData(const mxArray *prhs[]) {
    /**
     * @brief setupUserData extracts information from the matlab call and returns the corresponding UserData struct
     * @param[in] prhs: pointer to the array of input arguments @type mxArray
     * @return udata: struct containing all provided user data @type UserData
     */
    
    UserData udata; /* returned udata struct */
    realtype *plistdata; /* input for plist */
    realtype stldetdata; /* input for stldet */
    
    int ip;
    
    /* User udata structure */
    udata = (UserData) mxMalloc(sizeof *udata);
    if (udata == NULL) return(NULL);
    
    init_modeldims(udata);
    
    /* time */
    
    if (!prhs[0]) {
        mexErrMsgIdAndTxt("AMICI:mex:tout","No time vector provided!");
    }
    ts = mxGetPr(prhs[0]);
    
    nt = (int) mxGetM(prhs[0]) * mxGetN(prhs[0]);
    
    /* parameters */
    
    if (!prhs[1]) {
        mexErrMsgIdAndTxt("AMICI:mex:theta","No parameter vector provided!");
    }
    p = mxGetPr(prhs[1]);
    
    /* constants */
    
    if (!prhs[2]) {
        mexErrMsgIdAndTxt("AMICI:mex:kappa","No constant vector provided!");
    }
    k = mxGetPr(prhs[2]);
    
    if (!prhs[3]) {
        mexErrMsgIdAndTxt("AMICI:mex:options","No options provided!");
    }
    
    np = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);
    
    /* plist */
    if (!prhs[4]) {
        mexErrMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
    }
    
    if(prhs[4]) {
        plistdata = mxGetPr(prhs[4]);
    }
    
    plist = mxMalloc(np*sizeof(int));
    for (ip=0; ip<np; ip++) {
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
        idlist = mxGetData(mxGetProperty(prhs[3],0,"id")); \
    } else { \
        mexWarnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
        return(NULL); \
    }
    
    readOptionData(z2event)
    readOptionData(qpositivex)
    readOptionScalar(sensi,int)
    readOptionScalar(ism,int)
    readOptionScalar(sensi_meth,int)
    readOptionScalar(ordering,int)
    
    if(mxGetProperty(prhs[3], 0 ,"sx0")) { sx0data = mxGetPr(mxGetProperty(prhs[3], 0 ,"sx0"));} else { }
    if ((mxGetM(mxGetProperty(prhs[3], 0 ,"sx0")) * mxGetN(mxGetProperty(prhs[3], 0 ,"sx0")))>0) {
        /* check dimensions */
        if(mxGetN(mxGetProperty(prhs[3], 0 ,"sx0")) != np) { mexErrMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); }
        if(mxGetM(mxGetProperty(prhs[3], 0 ,"sx0")) != nx) { mexErrMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); }
        b_sx0 = TRUE;
    } else {
        b_sx0 = FALSE;
    }
    
    
    
    /* pbar */
    if (!prhs[5]) {
        mexErrMsgIdAndTxt("AMICI:mex:pbar","No parameter scales provided!");
    }
    
    pbar = mxGetPr(prhs[5]);
    
    /* xscale */
    if (!prhs[6]) {
        mexErrMsgIdAndTxt("AMICI:mex:xscale","No state scales provided!");
    }
    
    xbar = mxGetPr(prhs[6]);
    
    if (nx>0) {
        /* initialise temporary jacobian storage */
        tmp_J = NewSparseMat(nx,nx,nnz);
        M_tmp = mxMalloc(nx*nx*sizeof(realtype));
        dfdx_tmp = mxMalloc(nx*nx*sizeof(realtype));
    }
    if (sensi>0) {
        /* initialise temporary jacobian storage */
        tmp_dxdotdp = mxMalloc(nx*np*sizeof(realtype));
    }
    
    if (nw>0) {
        w_tmp = mxMalloc(nw*sizeof(realtype));
        dwdx_tmp = mxMalloc(ndwdx*sizeof(realtype));
        dwdp_tmp = mxMalloc(ndwdp*sizeof(realtype));
    }
    
    udata->am_nan_dxdotdp = FALSE;
    udata->am_nan_J = FALSE;
    udata->am_nan_JSparse = FALSE;
    udata->am_nan_xdot = FALSE;
    udata->am_nan_xBdot = FALSE;
    udata->am_nan_qBdot = FALSE;
    
    return(udata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

ReturnData setupReturnData(mxArray *plhs[], void *user_data, double *pstatus) {
    /**
     * setupReturnData initialises the return data struct
     * @param[in] plhs user input @type mxArray
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] pstatus pointer to the flag indicating the execution status @type double
     * @return rdata: return data struct @type ReturnData
     */
    ReturnData rdata; /* returned rdata struct */
    UserData udata; /** user udata */
    
    mxArray *mxsol;

    const char *field_names_sol[] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};
    mxArray *mxstatus;
    mxArray *mxllh;
    mxArray *mxsllh;
    mxArray *mxs2llh;
    mxArray *mxchi2;
    mxArray *mxt;
    mxArray *mxnumsteps;
    mxArray *mxnumrhsevals;
    mxArray *mxorder;
    mxArray *mxnumstepsS;
    mxArray *mxnumrhsevalsS;
    mxArray *mxrz;
    mxArray *mxz;
    mxArray *mxx;
    mxArray *mxy;
    mxArray *mxsrz;
    mxArray *mxsz;
    mxArray *mxsx;
    mxArray *mxsy;
    mxArray *mxs2rz;
    mxArray *mxsigmay;
    mxArray *mxssigmay;
    mxArray *mxsigmaz;
    mxArray *mxssigmaz;
    mxArray *mxxdot;
    mxArray *mxJ;
    mxArray *mxdydp;
    mxArray *mxdydx;
    mxArray *mxdxdotdp;

    mxArray *mxts;

    mwSize dimssx[] = {0,0,0};
    mwSize dimssy[] = {0,0,0};
    mwSize dimssz[] = {0,0,0};
    mwSize dimssrz[] = {0,0,0};
    mwSize dimss2rz[] = {0,0,0,0};
    mwSize dimssigmay[] = {0,0,0};
    mwSize dimssigmaz[] = {0,0,0};
    mwSize dimsssigmay[] = {0,0,0};
    mwSize dimsssigmaz[] = {0,0,0};
    
    /* this casting is necessary to ensure availability of accessor macros */
    udata = (UserData) user_data;
    
    /* Return rdata structure */
    rdata = (ReturnData) mxMalloc(sizeof *rdata);
    if (rdata == NULL) return(NULL);
    
    mxsol = mxCreateStructMatrix(1,1,29,field_names_sol);
    
    plhs[0] = mxsol;
    
    
    mxstatus = mxCreateDoubleMatrix(1,1,mxREAL);
    
    mxSetPr(mxstatus,pstatus);
    mxSetField(mxsol,0,"status",mxstatus);
    
    initField2(llh,1,1);
    initField2(chi2,1,1);
     
    mxts = mxCreateDoubleMatrix(nt,1,mxREAL);
    tsdata = mxGetPr(mxts);
    mxSetField(mxsol,0,"t",mxts);
    
    initField2(numsteps,nt,1);
    initField2(numrhsevals,nt,1);
    initField2(order,nt,1);
    if(sensi>0){
        initField2(numstepsS,nt,1);
        initField2(numrhsevalsS,nt,1);
    }
    if((nz>0) & (ne>0)){
        initField2(z,nmaxevent,nz);
        initField2(rz,nmaxevent,nztrue);
        initField2(sigmaz,nmaxevent,nz);
    }
    if(nx>0) {
        initField2(x,nt,nx);
        initField2(xdot,1,nx);
        initField2(J,nx,nx);
    }
    if(ny>0) {
        initField2(y,nt,ny);
        initField2(sigmay,nt,ny);
        if (sensi_meth == AMI_SS) {
            initField2(dydp,ny,np);
            initField2(dydx,ny,nx);
            initField2(dxdotdp,nx,np);
        }
    }
    if(sensi>0) {
        initField2(sllh,np,1);
        if (sensi_meth == AMI_FSA) {
            initField3(sx,nt,nx,np);
            if(ny>0) {
                initField3(sy,nt,ny,np);
                initField3(ssigmay,nt,ny,np);
            }
            if((nz>0) & (ne>0)){
                initField3(srz,nmaxevent,nztrue,np);
                if(sensi>1){
                    initField4(s2rz,nmaxevent,nztrue,np,np);
                }
                initField3(sz,nmaxevent,nz,np);
                initField3(ssigmaz,nmaxevent,nz,np);
            }
        }
        if (sensi_meth == AMI_ASA) {
            if(ny>0) {
                initField3(ssigmay,nt,ny,np);
            }
            if((nz>0) & (ne>0)){
                initField3(ssigmaz,nmaxevent,nz,np);
            }
        }
        if(sensi>1) {
            initField2(s2llh,np,np);
        }
    }
    
    return(rdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

ExpData setupExpData(const mxArray *prhs[], void *user_data) {
    /**
     * setupExpData initialises the experimental data struct
     * @param[in] prhs user input @type *mxArray
     * @param[in] user_data pointer to the user data struct @type UserData
     * @return edata: experimental data struct @type ExpData
     */
    
    int nmyt, nmyy, nysigmat, nysigmay; /* integers with problem dimensionality */
    int nmzt, nmzy, nzsigmat, nzsigmay; /* integers with problem dimensionality */
    
    char *errmsg;
    
    ExpData edata; /* returned rdata struct */
    UserData udata; /** user udata */
    
    errmsg = (char *)mxMalloc(200*sizeof(char));
    
    udata = (UserData) user_data;
    
    /* Return rdata structure */
    edata = (ExpData) mxMalloc(sizeof *edata);
    if (edata == NULL) return(NULL);
    
    if (!prhs[7]) {
        mexErrMsgIdAndTxt("AMICI:mex:data","No data provided!");
    }
    if (mxGetProperty(prhs[7], 0 ,"Y")) {
        my = mxGetPr(mxGetProperty(prhs[7], 0 ,"Y"));
        nmyy = (int) mxGetN(mxGetProperty(prhs[7], 0 ,"Y"));
        nmyt = (int) mxGetM(mxGetProperty(prhs[7], 0 ,"Y"));
    } else {
        mexErrMsgIdAndTxt("AMICI:mex:data:Y","Field Y not specified as field in data struct!");
    }
    
    if (mxGetProperty(prhs[7], 0 ,"Sigma_Y")) {
        ysigma = mxGetPr(mxGetProperty(prhs[7], 0 ,"Sigma_Y"));
        nysigmay = (int) mxGetN(mxGetProperty(prhs[7], 0 ,"Sigma_Y"));
        nysigmat = (int) mxGetM(mxGetProperty(prhs[7], 0 ,"Sigma_Y"));
    } else {
        mexErrMsgIdAndTxt("AMICI:mex:data:Sigma_Y","Field Sigma_Y not specified as field in data struct!");
    }
    if (mxGetProperty(prhs[7], 0 ,"Z")) {
        mz = mxGetPr(mxGetProperty(prhs[7], 0 ,"Z"));
        nmzy = (int) mxGetN(mxGetProperty(prhs[7], 0 ,"Z"));
        nmzt = (int) mxGetM(mxGetProperty(prhs[7], 0 ,"Z"));
    } else {
        mexErrMsgIdAndTxt("AMICI:mex:data:Z","Field Z not specified as field in data struct!");
    }
    
    if (mxGetProperty(prhs[7], 0 ,"Sigma_Z")) {
        zsigma = mxGetPr(mxGetProperty(prhs[7], 0 ,"Sigma_Z"));
        nzsigmay = (int) mxGetN(mxGetProperty(prhs[7], 0 ,"Sigma_Z"));
        nzsigmat = (int) mxGetM(mxGetProperty(prhs[7], 0 ,"Sigma_Z"));
    } else {
        mexErrMsgIdAndTxt("AMICI:mex:data:Sigma_Z","Field Sigma_Z not specified as field in data struct!");
    }
    
    if (nmyt != nt) {
        sprintf(errmsg,"Number of time-points in data matrix does (%i) not match provided time vector (%i)",nmyt,nt);
        mexErrMsgIdAndTxt("AMICI:mex:data:nty",errmsg);
    }
    
    if (nysigmat != nt) {
        sprintf(errmsg,"Number of time-points in data-sigma matrix (%i) does not match provided time vector (%i)",nysigmat,nt);
        mexErrMsgIdAndTxt("AMICI:mex:data:ntsdy",errmsg);
    }
    
    if (nmyy != nytrue) {
        sprintf(errmsg,"Number of observables in data matrix (%i) does not match model ny (%i)",nmyy,nytrue);
        mexErrMsgIdAndTxt("AMICI:mex:data:nyy",errmsg);
    }
    
    if (nysigmay != nytrue) {
        sprintf(errmsg,"Number of observables in data-sigma matrix (%i) does not match model ny (%i)",nysigmay,nytrue);
        mexErrMsgIdAndTxt("AMICI:mex:data:nysdy",errmsg);
    }
    
    if (nmzt != nmaxevent) {
        sprintf(errmsg,"Number of time-points in event matrix (%i) does not match provided nmaxevent (%i)",nmzt,nmaxevent);
        mexErrMsgIdAndTxt("AMICI:mex:data:nmaxeventnz",errmsg);
    }
    
    if (nzsigmat != nmaxevent) {
        sprintf(errmsg,"Number of time-points in event-sigma matrix (%i) does not match provided nmaxevent (%i)",nzsigmat,nmaxevent);
        mexErrMsgIdAndTxt("AMICI:mex:data:nmaxeventnsdz",errmsg);
    }
    
    if (nmzy != nztrue) {
        sprintf(errmsg,"Number of events in event matrix (%i) does not match provided nz (%i)",nmzy,nztrue);
        mexErrMsgIdAndTxt("AMICI:mex:data:nenz",errmsg);
    }
    
    if (nzsigmay != nztrue) {
        sprintf(errmsg,"Number of events in event-sigma matrix (%i) does not match provided nz (%i)",nzsigmay,nztrue);
        mexErrMsgIdAndTxt("AMICI:mex:data:nensdz",errmsg);
    }
    
    return(edata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void *setupAMI(int *status, void *user_data, void *temp_data) {
    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[in] temp_data pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block
     */
    void *ami_mem; /* pointer to ami memory block */
    bool error_corr = TRUE;
    int ip;
    int ix;
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    TempData tdata; /* user udata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;

    
    t = tstart;
    
    r = 0;
    g = 0;
    
    if (nx > 0) {
        
        /* allocate temporary objects */
        x = N_VNew_Serial(nx);
        x_old = N_VNew_Serial(nx);
        dx = N_VNew_Serial(nx); /* only needed for idas */
        dx_old = N_VNew_Serial(nx); /* only needed for idas */
        xdot = N_VNew_Serial(nx);
        xdot_old = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);
        
        if(ne>0) rootsfound = mxMalloc(ne*sizeof(int));
        if(ne>0) rootvals= mxMalloc(ne*sizeof(realtype));
        if(ne>0) rootidx = mxMalloc(nmaxevent*ne*ne*sizeof(int));
        if(ne>0) nroots = mxMalloc(ne*sizeof(int));
        if(ne>0) memset(nroots,0,ne*sizeof(int));
        if(ne>0) discs = mxMalloc(nmaxevent*ne*sizeof(realtype));
        if(ne>0) h = mxMalloc(ne*sizeof(realtype));
        if(ne>0) h_tmp = mxMalloc(ne*sizeof(realtype));
        
        if(ne>0) deltax = mxMalloc(nx*sizeof(realtype));
        if(ne>0) deltasx = mxMalloc(nx*np*sizeof(realtype));
        if(ne>0) deltaxB = mxMalloc(nx*sizeof(realtype));
        if(ne>0) deltaqB = mxMalloc(np*sizeof(realtype));
        
        if(ny>0) sigma_y = mxMalloc(ny*sizeof(realtype));
        if(ny>0) memset(sigma_y,0,ny*sizeof(realtype));
        if(ne>0) sigma_z = mxMalloc(nz*sizeof(realtype));
        if(ne>0) memset(sigma_z,0,nz*sizeof(realtype));
        
        
        /* initialise states */
        
        if (x == NULL) return(NULL);
        *status = fx0(x, udata);
        if (*status != AMI_SUCCESS) return(NULL);
        *status = fdx0(x, dx, udata); /* only needed for idas */
        if (*status != AMI_SUCCESS) return(NULL);
        
        /* initialise heaviside variables */
        
        initHeaviside(status,user_data,temp_data);
        
    }
    
    /* Create AMIS object */
    if (lmm>2||lmm<1) {
        mexErrMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        mexErrMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
    }
    ami_mem = AMICreate(lmm, iter);
    if (ami_mem == NULL) return(NULL);
    
    /* Initialize AMIS solver*/
    *status = wrap_init(ami_mem, x, dx, tstart);
    if (*status != AMI_SUCCESS) return(NULL);
    
    /* Specify integration tolerances */
    *status = AMISStolerances(ami_mem, RCONST(rtol), RCONST(atol));
    if(*status != AMI_SUCCESS) return(NULL);
    
    /* Set optional inputs */
    *status = AMISetErrHandlerFn(ami_mem);
    if(*status != AMI_SUCCESS) return(NULL);
    
    /* attaches userdata*/
    *status = AMISetUserData(ami_mem, udata);
    if(*status != AMI_SUCCESS) return(NULL);
    
    /* specify maximal number of steps */
    *status = AMISetMaxNumSteps(ami_mem, maxsteps);
    if(*status != AMI_SUCCESS) return(NULL);
    
    /* activates stability limit detection */
    *status = AMISetStabLimDet(ami_mem, stldet);
    if(*status != AMI_SUCCESS) return(NULL);
    
    if (ne>0) {
        /* activates root detection */
        *status = wrap_RootInit(ami_mem, udata);
        if(*status != AMI_SUCCESS) return(NULL);
    }
    
    /* Attach linear solver module */
    switch (linsol) {
            
            /* DIRECT SOLVERS */
            
        case AMI_DENSE:
            *status = AMIDense(ami_mem, nx);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetDenseJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case AMI_BAND:
            *status = AMIBand(ami_mem, nx, ubw, lbw);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetBandJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case AMI_LAPACKDENSE:
            mexErrMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* *status = CVLapackDense(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetDenseJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;
             
             break;*/
            
        case AMI_LAPACKBAND:
            
            mexErrMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* *status = CVLapackBand(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetBandJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;
             
             break;*/
            
        case AMI_DIAG:
            *status = AMIDiag(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMI_SPGMR:
            *status = AMISpgmr(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);

            break;
            
        case AMI_SPBCG:
            *status = AMISpbcg(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case AMI_SPTFQMR:
            *status = AMISptfqmr(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
            /* SPARSE SOLVERS */
            
        case AMI_KLU:
            *status = AMIKLU(ami_mem, nx, nnz);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetSparseJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
        
            *status = AMIKLUSetOrdering(ami_mem, ordering);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        default:
            mexErrMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");
            break;
    }
    
    if ( sensi >= 1) {
        
        dydx = mxMalloc(ny*nx*sizeof(realtype));
        memset(dydx,0,ny*nx*sizeof(realtype));
        dydp = mxMalloc(ny*np*sizeof(realtype));
        memset(dydp,0,ny*np*sizeof(realtype));
        
        dsigma_ydp = mxMalloc(ny*np*sizeof(realtype));
        memset(dsigma_ydp,0,ny*np*sizeof(realtype));
        if(ne>0) dsigma_zdp = mxMalloc(nz*np*sizeof(realtype));
        if(ne>0) memset(dsigma_zdp,0,nz*np*sizeof(realtype));
        
        if (sensi_meth == AMI_FSA) {
            
            if(nx>0) {
                
                /* allocate some more temporary storage */
                
                NVsx = N_VCloneVectorArray_Serial(np, x);
                sdx = N_VCloneVectorArray_Serial(np, x);
                if (NVsx == NULL) return(NULL);
                if (sdx == NULL) return(NULL);
                
                /* initialise sensitivities, this can either be user provided or come from the model definition */
                
                if(!b_sx0) {
                    *status = fsx0(NVsx, x, dx, udata);
                    if (*status != AMI_SUCCESS) return(NULL);
                } else {
                    for (ip=0; ip<np; ip++) {
                        sx_tmp = NV_DATA_S(NVsx[plist[ip]]);
                        for (ix=0; ix<nx; ix++) {
                            sx_tmp[ix] = sx0data[ix + nx*plist[ip]];
                        }
                    }
                }
                *status = fsdx0(sdx, x, dx, udata);
                if (*status != AMI_SUCCESS) return(NULL);
                
                /* Activate sensitivity calculations */
                
                *status = wrap_SensInit1(ami_mem, NVsx, sdx, udata);
                if (*status != AMI_SUCCESS) return(NULL);
                
                /* Set sensitivity analysis optional inputs */
                *status = AMISetSensParams(ami_mem, p, pbar, plist);
                if (*status != AMI_SUCCESS) return(NULL);
                
                *status = AMISetSensErrCon(ami_mem, error_corr);
                if (*status != AMI_SUCCESS) return(NULL);
                
                *status = AMISensEEtolerances(ami_mem);
                if (*status != AMI_SUCCESS) return(NULL);
            }
        }
        
        if (sensi_meth == AMI_ASA) {
            
            if(nx>0) {
                /* Allocate space for the adjoint computation */
                
                which = 0;
                
                if(ne>0) x_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);
                if(ne>0) xdot_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);
                if(ne>0) xdot_old_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);
                
                
                /* we always want N_d to be equal to the number of maximal steps, this prevents additional forward passes
                 and thus ensures correctness for systems with discontinuous right hand sides */
                
                *status = AMIAdjInit(ami_mem, maxsteps, interpType);
                if (*status != AMI_SUCCESS) return(NULL);
                
                llhS0 = mxMalloc(np*sizeof(realtype));
                memset(llhS0,0,np*sizeof(realtype));
                dgdp = mxMalloc(np*sizeof(realtype));
                memset(dgdp,0,np*sizeof(realtype));
                dgdx = mxMalloc(nx*nt*sizeof(realtype));
                memset(dgdx,0,nx*nt*sizeof(realtype));
                if (ne > 0) {
                    dzdp = mxMalloc(nz*np*sizeof(realtype));
                    memset(dzdp,0,nz*np*sizeof(realtype));
                    dzdx = mxMalloc(nz*nx*sizeof(realtype));
                    memset(dzdx,0,nz*nx*sizeof(realtype));
                }
                drdp = mxMalloc(np*nz*nmaxevent*sizeof(realtype));
                memset(drdp,0,np*nz*nmaxevent*sizeof(realtype));
                drdx = mxMalloc(nx*nz*nmaxevent*sizeof(realtype));
                memset(drdx,0,nx*nz*nmaxevent*sizeof(realtype));
            }
        }
        
        
        
    }
    
    id = N_VNew_Serial(nx);
    id_tmp = NV_DATA_S(id);
    memcpy(id_tmp,idlist,nx*sizeof(realtype));
    
    *status = AMISetId(ami_mem, id);
    if (*status != AMI_SUCCESS) return(NULL);
    
    *status = AMISetSuppressAlg(ami_mem, TRUE);
    if (*status != AMI_SUCCESS) return(NULL);

    
    return(ami_mem);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void setupAMIB(int *status,void *ami_mem, void *user_data, void *temp_data) {
    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[in] temp_data pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block for the backward problem
     */
    /* this casting is necessary to ensure availability of accessor macros */
    int ix;
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    

    
    xB = N_VNew_Serial(nx);
    xB_old = N_VNew_Serial(nx);
    
    dxB = N_VNew_Serial(nx);
    
    xQB = N_VNew_Serial(np);
    xQB_old = N_VNew_Serial(np);
    
    /* write initial conditions */
    if (xB == NULL) return;
    xB_tmp = NV_DATA_S(xB);
    memset(xB_tmp,0,sizeof(realtype)*nx);
    for (ix=0; ix<nx; ix++) {
        xB_tmp[ix] += dgdx[nt-1+ix*nt];
    }
    
    if (dxB == NULL) return;
    dxB_tmp = NV_DATA_S(dxB);
    memset(dxB_tmp,0,sizeof(realtype)*nx);
    
    if (xQB == NULL) return;
    xQB_tmp = NV_DATA_S(xQB);
    memset(xQB_tmp,0,sizeof(realtype)*np);
    
    /* create backward problem */
    if (lmm>2||lmm<1) {
        mexErrMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        mexErrMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
    }
    /* allocate memory for the backward problem */
    *status = AMICreateB(ami_mem, lmm, iter, &which);
    if (*status != AMI_SUCCESS) return;
    
    
    /* initialise states */
    *status = wrap_binit(ami_mem, which, xB, dxB, t);
    if (*status != AMI_SUCCESS) return;
    
    /* specify integration tolerances for backward problem */
    *status = AMISStolerancesB(ami_mem, which, RCONST(rtol), RCONST(atol));
    if(*status != AMI_SUCCESS) return;
    
    /* Attach user data */
    *status = AMISetUserDataB(ami_mem, which, udata);
    if(*status != AMI_SUCCESS) return;
    
    /* Number of maximal internal steps */
    *status = AMISetMaxNumStepsB(ami_mem, which, 100*maxsteps);
    if(*status != AMI_SUCCESS) return;
    
    switch (linsol) {
            
            /* DIRECT SOLVERS */
            
        case AMI_DENSE:
            *status = AMIDenseB(ami_mem, which, nx);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetDenseJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case AMI_BAND:
            *status = AMIBandB(ami_mem, which, nx, ubw, lbw);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetBandJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case AMI_LAPACKDENSE:
            
            /* #if SUNDIALS_BLAS_LAPACK
             *status = CVLapackDenseB(ami_mem, which, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetDenseJacFnB(ami_mem, which);
             if (*status != AMI_SUCCESS) return;
             #else*/
            mexErrMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMI_LAPACKBAND:
            
            
            /* #if SUNDIALS_BLAS_LAPACK
             *status = CVLapackBandB(ami_mem, which, nx, ubw, lbw);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetBandJacFnB(ami_mem, which);
             if (*status != AMI_SUCCESS) return;
             #else*/
            mexErrMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
            break;
            break;
            
        case AMI_DIAG:
            *status = AMIDiagB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetDenseJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case AMI_SPGMR:
            *status = AMISpgmrB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case AMI_SPBCG:
            *status = AMISpbcgB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case AMI_SPTFQMR:
            *status = AMISptfqmrB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
            /* SPARSE SOLVERS */
            
        case AMI_KLU:
            *status = AMIKLUB(ami_mem, which, nx, nnz);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetSparseJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            *status = AMIKLUSetOrderingB(ami_mem, which, ordering);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        default:
            break;
    }
    
    /* Initialise quadrature calculation */
    *status = wrap_qbinit(ami_mem, which, xQB);
    if (*status != AMI_SUCCESS) return;
    
    /* Enable Quadrature Error Control */
    *status = AMISetQuadErrConB(ami_mem, which, TRUE);
    if (*status != AMI_SUCCESS) return;
    
    *status = AMIQuadSStolerancesB(ami_mem, which, RCONST(rtol), RCONST(atol));
    if(*status != AMI_SUCCESS) return;
    
    *status = AMISetStabLimDetB(ami_mem, which, stldet); /* activates stability limit detection */
    if(*status != AMI_SUCCESS) return;
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDataSensisFSA(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ip;
    int iy;
    int ix;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    for(ip=0; ip < np; ip++) {
        if(nx>0) {
            if(ts[it] > tstart) {
                *status = AMIGetSens(ami_mem, &t, NVsx);
                if (*status != AMI_SUCCESS) return;
            }
            
            sx_tmp = NV_DATA_S(NVsx[ip]);
            for(ix=0; ix < nx; ix++) {
                sxdata[(ip*nx + ix)*nt + it] = sx_tmp[ix];
            }
        }
    }
    for (iy=0; iy<ny; iy++) {
        if (mxIsNaN(ysigma[iy*nt+it])) {
            *status = fdsigma_ydp(t,dsigma_ydp,udata);
            if (*status != AMI_SUCCESS) return;
        } else {
            for (ip=0; ip<np; ip++) {
                dsigma_ydp[ip*ny+iy] = 0;
            }
        }
        for (ip=0; ip<np; ip++) {
            ssigmaydata[it + nt*(ip*ny+iy)] = dsigma_ydp[ip*ny+iy];
        }
    }
    fdydx(ts[it],it,dydx,x,udata);
    fdydp(ts[it],it,dydp,x,udata);
    fsy(ts[it],it,sydata,dydx,dydp,NVsx,udata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDataSensisASA(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getDataSensisASA extracts data information for adjoint sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int iy;
    int ip;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    *status = fdydx(ts[it],it,dydx,x,udata);
    if (*status != AMI_SUCCESS) return;
    *status = fdydp(ts[it],it,dydp,x,udata);
    if (*status != AMI_SUCCESS) return;
    for (iy=0; iy<ny; iy++) {
        if (mxIsNaN(ysigma[iy*nt+it])) {
            *status = fdsigma_ydp(t,dsigma_ydp,udata);
            if (*status != AMI_SUCCESS) return;
        } else {
            for (ip=0; ip<np; ip++) {
                dsigma_ydp[ip*ny+iy] = 0;
            }
        }
        for (ip=0; ip<np; ip++) {
            ssigmaydata[it + nt*(ip*ny+iy)] = dsigma_ydp[ip*ny+iy];
        }
    }
    fdJydp(ts[it],it,dgdp,ydata,x,dydp,my,sigma_y,dsigma_ydp,udata);
    fdJydx(ts[it],it,dgdx,ydata,x,dydx,my,sigma_y,udata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDataOutput(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int iy;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    *status = fy(ts[it],it,ydata,x,udata);
    if (*status != AMI_SUCCESS) return;
    
    for (iy=0; iy<ny; iy++) {
        /* extract the value for the standard deviation, if the data value is NaN, use
         the parameter value. Store this value in the return struct */
        if (mxIsNaN(ysigma[iy*nt+it])) {
            *status =fsigma_y(t,sigma_y,udata);
            if (*status != AMI_SUCCESS) return;
            
        } else {
            sigma_y[iy] = ysigma[iy*nt+it];
        }
        sigmaydata[iy*nt+it] = sigma_y[iy];
    }
    fJy(t,it,&g,ydata,x,my,sigma_y,udata);
    if (sensi >= 1) {
        if (sensi_meth == AMI_FSA) {
            getDataSensisFSA(status, it, ami_mem, udata, rdata, edata, tdata);
        }
        if (sensi_meth == AMI_ASA) {
            getDataSensisASA(status, it, ami_mem, udata, rdata, edata, tdata);
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisFSA(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    TempData tdata; /* temp data */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    tdata = (TempData) temp_data;

    *status = fsz(t,ie,nroots,szdata,x,NVsx,udata);
    if (*status != AMI_SUCCESS) return;

}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisFSA_tf(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getEventSensisFSA_tf extracts event information for forward sensitivity
     *     analysis for events that happen at the end of the considered interval
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    TempData tdata; /* temp data */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    tdata = (TempData) temp_data;
    
    *status = fsz_tf(t,ie,nroots,szdata,x,NVsx,udata);
    if (*status != AMI_SUCCESS) return;
    
    *status = fsroot(t,ie,nroots,srzdata,x,NVsx,udata);
    if (*status != AMI_SUCCESS) return;
    
    if(sensi>1) {
        *status = fs2root(t,ie,nroots,s2rzdata,x,NVsx,udata);
        if (*status != AMI_SUCCESS) return;
    }
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisASA(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventSensisASA extracts event information for adjoint sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int ip;
    int ix;
    int iz;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    for (iz=0; iz<nztrue; iz++) {
        if( z2event[iz] == ie ){
            if(!mxIsNaN(mz[iz*nmaxevent+nroots[ie]])) {
                *status = fdzdp(t,ie,dzdp,x,udata);
                if (*status != AMI_SUCCESS) return;
                *status = fdzdx(t,ie,dzdx,x,udata);
                if (*status != AMI_SUCCESS) return;
                /* extract the value for the standard deviation, if the data value is NaN, use
                 the parameter value. Store this value in the return struct */
                if (mxIsNaN(zsigma[nroots[ie] + nmaxevent*iz])) {
                    *status = fsigma_z(t,ie,sigma_z,udata);
                    if (*status != AMI_SUCCESS) return;
                    *status = fdsigma_zdp(t,ie,dsigma_zdp,udata);
                    if (*status != AMI_SUCCESS) return;
                } else {
                    for (ip=0; ip<np; ip++) {
                        dsigma_zdp[iz+nz*ip] = 0;
                    }
                    sigma_z[iz] = zsigma[nroots[ie] + nmaxevent*iz];
                }
                sigmazdata[nroots[ie] + nmaxevent*iz] = sigma_z[iz];
                for (ip=0; ip<np; ip++) {
                    ssigmazdata[nroots[ie] + nmaxevent*(iz+nz*ip)] = dsigma_zdp[iz+nz*ip];
                }
                
                for (ip=0; ip<np; ip++) {
                    drdp[nroots[ie] + nmaxevent*ip] += dsigma_zdp[ip*ne+ie]/sigma_z[iz] + ( dzdp[ip +np*iz]* ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] ) )/pow( sigma_z[ie] , 2) - dsigma_zdp[ip*ne+ie]*pow( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] ,2)/pow( sigma_z[iz] , 3);
                }
                for (ix=0; ix<nx; ix++) {
                    drdx[nroots[ie] + nmaxevent*ix] += ( dzdx[ip + nx*iz] * ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] ) )/pow( sigma_z[iz] , 2);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSigma(int *status, int ie, int iz, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventSigma extracts fills sigma_z either from the user defined function or from user input
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie event type index @type int
     * @param[in] iz event output index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    /* extract the value for the standard deviation, if the data value is NaN, use
     the parameter value. Store this value in the return struct */
    if (mxIsNaN(zsigma[nroots[ie] + nmaxevent*iz])) {
        *status = fsigma_z(t,ie,sigma_z,udata);
        if (*status != AMI_SUCCESS) return;
    } else {
        sigma_z[iz] = zsigma[nroots[ie] + nmaxevent*iz];
    }
    sigmazdata[nroots[ie] + nmaxevent*iz] = sigma_z[iz];
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventObjective(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventObjective updates the objective function on the occurence of an event
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie event type index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int iz,ip;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    for (iz=0; iz<nztrue; iz++) {
        if(z2event[iz] == ie) {
            getEventSigma(status, ie, iz, ami_mem, user_data, return_data, exp_data, temp_data);
            if(!mxIsNaN(mz[iz*nmaxevent+nroots[ie]])) {
                r += 0.5*log(2*pi*pow(zsigma[nroots[ie] + nmaxevent*iz],2)) + 0.5*pow( ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] )/zsigma[iz] , 2);
                *chi2data += pow( ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] )/zsigma[iz] , 2);
            }
        }
    }
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventOutput(int *status, realtype *tlastroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventOutput extracts output information for events
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] tlastroot timepoint of last occured event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int iz;
    int ie;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    
    /* EVENT OUTPUT */
    for (ie=0; ie<ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
        if (nroots[ie]<nmaxevent) {
            if(rootsfound[ie] != 0) {
                *status = fz(t,ie,nroots,zdata,x,udata);
                if (*status != AMI_SUCCESS) return;
                
                for (iz=0; iz<nztrue; iz++) {
                    if(z2event[iz] == ie) {
                        getEventSigma(status, ie, iz, ami_mem,user_data,return_data,exp_data,temp_data);
                        if (*status != AMI_SUCCESS) return;
                    }
                }

                getEventObjective(status, ie, ami_mem, user_data, return_data, exp_data, temp_data);
                if (*status != AMI_SUCCESS) return;                
                
                if (sensi >= 1) {
                    if(sensi_meth == AMI_ASA) {
                        getEventSensisASA(status, ie, ami_mem, udata, rdata, edata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    } else {
                        getEventSensisFSA(status, ie, ami_mem, udata, rdata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    }
                }
                
                nroots[ie]++;
            }
        }
    }
    return;
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void fillEventOutput(int *status, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * fillEventOutput fills missing roots at last timepoint
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ie;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    froot(t,x,dx,rootvals,user_data);
    
    
    /* EVENT OUTPUT */
    if (nztrue>0) {
        for (ie=0; ie<ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
            if (nroots[ie]<nmaxevent) {
                *status = fz(t,ie,nroots,zdata,x,udata);
                if (*status != AMI_SUCCESS) return;
                
                
                rzdata[nroots[ie]+ie] = rootvals[ie];
                
                
                getEventObjective(status, ie, ami_mem, user_data, return_data, exp_data, temp_data);
                if (*status != AMI_SUCCESS) return;
                
                if (sensi >= 1) {
                    if(sensi_meth == AMI_ASA) {
                        getEventSensisASA(status, ie, ami_mem, udata, rdata, edata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    } else {
                        getEventSensisFSA_tf(status, ie, ami_mem, udata, rdata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    }
                }
                
                nroots[ie]++;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleDataPoint(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ix;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    tsdata[it] = ts[it];
    x_tmp = NV_DATA_S(x);
    for (ix=0; ix<nx; ix++) {
        xdata[it+nt*ix] = x_tmp[ix];
    }
    
    if (it == nt-1) {
        if( sensi_meth == AMI_SS) {
            
            *status = fdxdotdp(t,dxdotdpdata,x,dx,udata);
            if (*status != AMI_SUCCESS) return;
            *status = fdydp(ts[it],it,dydpdata,x,udata);
            if (*status != AMI_SUCCESS) return;
            *status = fdydx(ts[it],it,dydxdata,x,udata);
            if (*status != AMI_SUCCESS) return;
        }
    }
    
    if(ts[it] > tstart) {
        getDiagnosis(status, it, ami_mem, udata, rdata);
        if (*status != AMI_SUCCESS) return;
    }
    
    getDataOutput(status, it, ami_mem, udata, rdata, edata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleDataPointB(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points for the backward problems
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ix;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    xB_tmp = NV_DATA_S(xB);
    for (ix=0; ix<nx; ix++) {
        xB_tmp[ix] += dgdx[it+ix*nt];
    }
    getDiagnosisB(status,it,ami_mem,user_data,return_data,temp_data);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleEvent(int *status, int *iroot, realtype *tlastroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data, int seflag) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] iroot index of event @type int
     * @param[out] tlastroot pointer to the timepoint of the last event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int ie;
    int secondevent = 0;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    /* store heaviside information at event occurence */
    froot(t,x,dx,rootvals,user_data);
    for (ie = 0; ie<ne; ie++) {
        h_tmp[ie] = rootvals[ie];
    }
    
    if (seflag == 0) {
        *status = AMIGetRootInfo(ami_mem, rootsfound);
        if (*status != AMI_SUCCESS) return;
    }
    
    if (*iroot<nmaxevent*ne) {
        for (ie=0; ie<ne; ie++) {
            rootidx[*iroot*ne + ie] = rootsfound[ie];
        }
    }
    
    /* only extract in the first event fired */
    if (seflag == 0) {
        if(sensi >= 1){
            if (sensi_meth == AMI_FSA) {
                *status = AMIGetSens(ami_mem, &t, NVsx);
                if (*status != AMI_SUCCESS) return;
            }
        }
    }
    
    /* only check this in the first event fired, otherwise this will always be true */
    if (seflag == 0) {
        if (t == *tlastroot) {
            /* we are stuck in a root => turn off rootfinding */
            /* at some point we should find a more intelligent solution here, and turn on rootfinding again after some time */
            AMIRootInit(ami_mem, 0, NULL);
        }
        *tlastroot = t;
    }
    
    getEventOutput(status, tlastroot, ami_mem, udata, rdata, edata, tdata);
    if (*status != AMI_SUCCESS) return;
    
    /* if we need to do forward sensitivities later on we need to store the old x and the old xdot */
    if(sensi >= 1){
        if (sensi_meth == AMI_FSA) {
            /* store x and xdot to compute jump in sensitivities */
            N_VScale(1.0,x,x_old);
            
            *status = fxdot(t,x,dx,xdot,udata);
            N_VScale(1.0,xdot,xdot_old);
            N_VScale(1.0,dx,dx_old);
        }
        
        if (sensi_meth == AMI_ASA) {
            /* store x to compute jump in discontinuity */
            N_VScale(1.0,x,x_disc[*iroot]);
            N_VScale(1.0,xdot,xdot_disc[*iroot]);
            N_VScale(1.0,xdot_old,xdot_old_disc[*iroot]);
        }
    }
    
    updateHeaviside(status, udata, tdata);
    if (*status != AMI_SUCCESS) return;
    
    applyEventBolus(status, ami_mem, udata, tdata);
    if (*status != AMI_SUCCESS) return;
    
    if (*iroot<nmaxevent*ne) {
        discs[*iroot] = t;
        (*iroot)++;
    } else {
        mexWarnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT","Event was recorded but not reported as the number of occured events exceeded (nmaxevents)*(number of events in model definition)!");
        *status = AMIReInit(ami_mem, t, x, dx); /* reinitialise so that we can continue in peace */
        return;
    }
    
    if(sensi >= 1){
        if (sensi_meth == AMI_FSA) {
            
            /* compute the new xdot  */
            *status = fxdot(t,x,dx,xdot,udata);
            if (*status != AMI_SUCCESS) return;
            
            applyEventSensiBolusFSA(status, ami_mem, udata, tdata);
            if (*status != AMI_SUCCESS) return;
        }
    }
    
    /* check whether we need to fire a secondary event */
    froot(t,x,dx,rootvals,user_data);
    for (ie = 0; ie<ne; ie++) {
        if(h_tmp[ie] != rootvals[ie]) {
            if (h_tmp[ie]<rootvals[ie]) {
                rootsfound[ie] = 1;
            } else {
                rootsfound[ie] = -1;
            }
            secondevent++;
        } else {
            rootsfound[ie] = 0;
        }
    }
    /* fire the secondary event */
    if(secondevent>0) {
        handleEvent(status, iroot, tlastroot, ami_mem, udata, rdata, edata, tdata, secondevent);
    }
    
    /* only reinitialise in the first event fired */
    if (seflag == 0) {
        *status = AMIReInit(ami_mem, t, x, dx);
        if (*status != AMI_SUCCESS) return;
        
        /* make time derivative consistent */
        *status = AMICalcIC(ami_mem, t);
        if (*status != AMI_SUCCESS) return;
    }
    
    if(sensi >= 1){
        if (sensi_meth == AMI_FSA) {
            if(sensi >= 1){
                if (sensi_meth == AMI_FSA) {
                    if (seflag == 0) {
                        *status = AMISensReInit(ami_mem, ism, NVsx, sdx);
                        if (*status != AMI_SUCCESS) return;
                    }
                }
            }
        }
    }
    
    return;
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleEventB(int *status, int iroot, void *ami_mem, void  *user_data, void *temp_data) {
    /**
     * handleEventB executes everything necessary for the handling of events for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] iroot index of event @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return cv_status updated status flag @type int
     */
    
    int ie;
    int ix;
    int ip;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    /* store current values */
    N_VScale(1.0,xB,xB_old);
    N_VScale(1.0,xQB,xQB_old);
    
    xB_tmp = NV_DATA_S(xB);
    xQB_tmp = NV_DATA_S(xQB);
    
    for (ie=0; ie<ne; ie++) {
        
        if (rootidx[iroot*ne + ie] != 0) {
            *status = fdeltaqB(t,ie,deltaqB,x_disc[iroot],xB_old,xQB_old,xdot_disc[iroot],xdot_old_disc[iroot],udata);
            if (*status != AMI_SUCCESS) return;
            *status = fdeltaxB(t,ie,deltaxB,x_disc[iroot],xB_old,xdot_disc[iroot],xdot_old_disc[iroot],udata);
            if (*status != AMI_SUCCESS) return;
            
            for (ix=0; ix<nx; ix++) {
                xB_tmp[ix] += deltaxB[ix];
                if (nz>0) {
                    xB_tmp[ix] += drdx[nroots[ie] + nmaxevent*ix];
                }
            }
            
            for (ip=0; ip<np; ip++) {
                xQB_tmp[ip] += deltaqB[ip];
            }
            nroots[ie]--;
            
        }
    }
    
    updateHeavisideB(status, iroot, udata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, void *user_data) {
    /**
     * getTnext computes the next timepoint to integrate to. This is the maximum of
     * tdata and troot but also takes into account if it<0 or iroot<0 where these expressions
     * do not necessarily make sense
     *
     * @param[in] troot timepoint of next event @type realtype
     * @param[in] iroot index of next event @type int
     * @param[in] tdata timepoint of next data point @type realtype
     * @param[in] it index of next data point @type int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @return tnext next timepoint @type realtype
     */
    
    realtype tnext;
    
    UserData udata; /* user udata */
    udata = (UserData) user_data;
    
    if (it<0) {
        tnext = troot[iroot];
    } else {
        if (iroot<0) {
            tnext = tdata[it];
        } else {
            if (ne>0) {
                if (troot[iroot]>tdata[it]) {
                    tnext = troot[iroot];
                } else {
                    tnext = tdata[it];
                }
            } else {
                tnext = tdata[it];
            }
        }
    }
    
    return(tnext);
    
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void applyEventBolus(int *status, void *ami_mem, void  *user_data, void *temp_data) {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ix;
    int ie;

    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    for (ie=0; ie<ne; ie++){
        if(rootsfound[ie] != 0) {
            *status = fdeltax(t,ie,deltax,x,xdot,xdot_old,user_data);
            
            x_tmp = NV_DATA_S(x);
            for (ix=0; ix<nx; ix++) {
                x_tmp[ix] += deltax[ix];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void applyEventSensiBolusFSA(int *status, void *ami_mem, void  *user_data, void *temp_data) {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current sensitivities
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ix;
    int ip;
    int ie;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    for (ie=0; ie<ne; ie++){
        if(rootsfound[ie] != 0) {
            *status = fdeltasx(t,ie,deltasx,x_old,xdot,xdot_old,NVsx,user_data);
            
            for (ip=0; ip<np; ip++) {
                sx_tmp = NV_DATA_S(NVsx[plist[ip]]);
                for (ix=0; ix<nx; ix++) {
                    sx_tmp[ix] += deltasx[ix + nx*ip];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void initHeaviside(int *status, void  *user_data, void *temp_data) {
    /**
     * initHeaviside initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ie;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    froot(t,x,dx,rootvals,user_data);
    
    for (ie = 0; ie<ne; ie++) {
        if (rootvals[ie]<0) {
            h[ie] = 0.0;
        } else {
            h[ie] = 1.0;
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void updateHeaviside(int *status, void  *user_data, void *temp_data) {
    /**
     * updateHeaviside updates the heaviside variables h on event occurences
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ie;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    /* rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */
    
    for (ie = 0; ie<ne; ie++) {
        h[ie] += rootsfound[ie];
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void updateHeavisideB(int *status, int iroot, void  *user_data, void *temp_data) {
    /**
     * updateHeavisideB updates the heaviside variables h on event occurences for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] iroot discontinuity occurance index @type int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ie;
    
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    /* rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */
    
    for (ie = 0; ie<ne; ie++) {
        h[ie] -= rootidx[iroot*ne + ie];
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDiagnosis(int *status,int it, void *ami_mem, void  *user_data, void *return_data) {
    /**
     * getDiagnosis extracts diagnosis information from solver memory block and writes them into the return data struct
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @return void
     */
    long int numsteps;
    long int numrhsevals;
    int order;

    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    
    *status = AMIGetNumSteps(ami_mem, &numsteps);
    if (*status != AMI_SUCCESS) return;
    numstepsdata[it] = (realtype)numsteps;
   
    *status = AMIGetNumRhsEvals(ami_mem, &numrhsevals);
    if (*status != AMI_SUCCESS) return;
    numrhsevalsdata[it] = (realtype)numrhsevals;
   
    *status = AMIGetLastOrder(ami_mem, &order);
    if (*status != AMI_SUCCESS) return;
    orderdata[it] = (realtype)order;

}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDiagnosisB(int *status,int it, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getDiagnosisB extracts diagnosis information from solver memory block and writes them into the return data struct for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    long int numsteps;
    long int numrhsevals;
    
    void *ami_memB;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    tdata = (TempData) temp_data;
    
    ami_memB = AMIGetAdjBmem(ami_mem, which);
    
    *status = AMIGetNumSteps(ami_memB, &numsteps);
    if (*status != AMI_SUCCESS) return;
    numstepsSdata[it] = (realtype)numsteps;
    
    *status = AMIGetNumRhsEvals(ami_memB, &numrhsevals);
    if (*status != AMI_SUCCESS) return;
    numrhsevalsSdata[it] = (realtype)numrhsevals;
    
}