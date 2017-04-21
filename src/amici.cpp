/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <math.h>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

#include "wrapfunctions.h" /* user functions */
#include <include/amici.h> /* amici functions */
#include <include/symbolic_functions.h>

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
#ifndef AMICI_WITHOUT_MATLAB
#define initField2(FIELD,D1,D2) \
mxArray *mx ## FIELD; \
mx ## FIELD = mxCreateDoubleMatrix(D1,D2,mxREAL); \
FIELD ## data = mxGetPr(mx ## FIELD); \
mxSetField(mxsol,0,#FIELD,mx ## FIELD)
#else
#define initField2(FIELD,D1,D2) \
double *mx ## FIELD; \
mx ## FIELD = new double[D1 * D2](); \
FIELD ## data = mx ## FIELD;
#endif

/**
 * @ brief initialise 3D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 */
#ifndef AMICI_WITHOUT_MATLAB
#define initField3(FIELD,D1,D2,D3) \
mxArray *mx ## FIELD; \
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
mx ## FIELD = mxCreateNumericArray(3,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
FIELD ## data = mxGetPr(mx ## FIELD); \
mxSetField(mxsol,0,#FIELD,mx ## FIELD)
#else
#define initField3(FIELD,D1,D2,D3) \
double *mx ## FIELD; \
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
mx ## FIELD = new double[D1 * D2 * D3](); \
FIELD ## data = mx ## FIELD;
#endif

/**
 * @ brief initialise 4D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 * @ param D4 number of elements in the fourth dimension of the tensor
 */
#ifndef AMICI_WITHOUT_MATLAB
#define initField4(FIELD,D1,D2,D3,D4) \
mxArray *mx ## FIELD; \
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
dims ## FIELD[3]=D4; \
mx ## FIELD = mxCreateNumericArray(4,dims ## FIELD,mxDOUBLE_CLASS,mxREAL); \
FIELD ## data = mxGetPr(mx ## FIELD); \
mxSetField(mxsol,0,#FIELD,mx ## FIELD)
#else
#define initField4(FIELD,D1,D2,D3,D4) \
double *mx ## FIELD; \
dims ## FIELD[0]=D1; \
dims ## FIELD[1]=D2; \
dims ## FIELD[2]=D3; \
dims ## FIELD[3]=D4; \
mx ## FIELD = new double[D1 * D2 * D3 * D4](); \
FIELD ## data = mx ## FIELD;
#endif

/**
 * @ brief extract information from a property of a matlab class (scalar)
 * @ param OPTION name of the property
 * @ param TYPE class to which the information should be cast
 */
#ifndef AMICI_WITHOUT_MATLAB
#define readOptionScalar(OPTION,TYPE) \
if(mxGetProperty(prhs[3],0,#OPTION)){ \
OPTION = (TYPE)mxGetScalar(mxGetProperty(prhs[3],0,#OPTION)); \
} else { \
warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
return(NULL); \
}
#endif

/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */
#ifndef AMICI_WITHOUT_MATLAB
#define readOptionData(OPTION) \
if(mxGetProperty(prhs[3],0,#OPTION)){ \
OPTION = (double *) mxGetData(mxGetProperty(prhs[3],0,#OPTION)); \
} else { \
warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
return(NULL); \
}
#endif

#ifdef AMICI_WITHOUT_MATLAB
    typedef double mxArray;
    void unscaleParameters(UserData *udata);
    void applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata);
#endif

/** return value for successful execution */
#define AMI_SUCCESS               0



#ifndef AMICI_WITHOUT_MATLAB
UserData *setupUserData(const mxArray *prhs[]) {
    /**
     * @brief setupUserData extracts information from the matlab call and returns the corresponding UserData struct
     * @param[in] prhs: pointer to the array of input arguments @type mxArray
     * @return udata: struct containing all provided user data @type UserData
     */

    UserData *udata; /* returned udata *struct */
    realtype *plistdata; /* input for plist */

    int ip;

    /* User udata structure */
    udata = new UserData();
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

    np = (int) mxGetM(prhs[4]) * mxGetN(prhs[4]);

    /* plist */
    if (!prhs[4]) {
        errMsgIdAndTxt("AMICI:mex:plist","No parameter list provided!");
    }

    if(prhs[4]) {
        plistdata = mxGetPr(prhs[4]);
    }

    plist = new int[np]();
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
        idlist = (double *) mxGetData(mxGetProperty(prhs[3],0,"id")); \
    } else { \
        warnMsgIdAndTxt("AMICI:mex:OPTION","Provided options are not of class amioption!"); \
        return(NULL); \
    }

    readOptionData(z2event)
    readOptionData(qpositivex)
    readOptionScalar(sensi,int)
    readOptionScalar(ism,int)
    readOptionScalar(sensi_meth,int)
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
    if (!prhs[7]) {
        b_x0 = FALSE;
        b_sx0 = FALSE;
    } else {
        if(mxGetField(prhs[7], 0 ,"x0")) {
            x0data = mxGetPr(mxGetField(prhs[7], 0 ,"x0"));
            if ((mxGetM(mxGetField(prhs[7], 0 ,"x0")) * mxGetN(mxGetField(prhs[7], 0 ,"x0")))>0) {
                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"x0")) != 1) { errMsgIdAndTxt("AMICI:mex:x0","Number of rows in x0 field must be equal to 1!"); }
                if(mxGetM(mxGetField(prhs[7], 0 ,"x0")) != nx) { errMsgIdAndTxt("AMICI:mex:x0","Number of columns in x0 field does not agree with number of model states!"); }
                b_x0 = TRUE;
            } else {
                b_x0 = FALSE;
            }
        } else {
            b_x0 = FALSE;
        }

        if(mxGetField(prhs[7], 0 ,"sx0")) {
            sx0data = mxGetPr(mxGetField(prhs[7], 0 ,"sx0"));
            if ((mxGetM(mxGetField(prhs[7], 0 ,"sx0")) * mxGetN(mxGetField(prhs[7], 0 ,"sx0")))>0) {
                /* check dimensions */
                if(mxGetN(mxGetField(prhs[7], 0 ,"sx0")) != np) { errMsgIdAndTxt("AMICI:mex:sx0","Number of rows in sx0 field does not agree with number of model parameters!"); }
                if(mxGetM(mxGetField(prhs[7], 0 ,"sx0")) != nx) { errMsgIdAndTxt("AMICI:mex:sx0","Number of columns in sx0 field does not agree with number of model states!"); }
                b_sx0 = TRUE;
            } else {
                b_sx0 = FALSE;
            }
        } else {
            b_sx0 = FALSE;
        }
    }


    if (nx>0) {
        /* initialise temporary jacobian storage */
        tmp_J = SparseNewMat(nx,nx,nnz,CSC_MAT);
        M_tmp = new realtype[nx*nx]();
        dfdx_tmp = new realtype[nx*nx]();
    }
    if (sensi>0) {
        /* initialise temporary dxdotdp storage */
        tmp_dxdotdp = new realtype[nx*np]();
    }
    if (ne>0) {
        /* initialise temporary stau storage */
        stau_tmp = new realtype[np]();
    }

    w_tmp = new realtype[nw]();
    dwdx_tmp = new realtype[ndwdx]();
    dwdp_tmp = new realtype[ndwdp]();

    return(udata);
}
#endif

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

#ifndef AMICI_WITHOUT_MATLAB
ReturnData *setupReturnData(mxArray *plhs[], UserData *udata, double *pstatus) {
    /**
     * setupReturnData initialises the return data struct
     * @param[in] plhs user input @type mxArray
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] pstatus pointer to the flag indicating the execution status @type double
     * @return rdata: return data struct @type ReturnData
     */
    ReturnData *rdata; /* returned rdata struct */

    mxArray *mxsol;

    const char *field_names_sol[] = {"status","llh","sllh","s2llh","chi2","t","numsteps","numrhsevals","order","numstepsS","numrhsevalsS","rz","z","x","y","srz","sz","sx","sy","s2rz","sigmay","ssigmay","sigmaz","ssigmaz","xdot","J","dydp","dydx","dxdotdp"};
    mxArray *mxstatus;

    mxArray *mxts;

    mwSize dimssx[] = {0,0,0};
    mwSize dimssy[] = {0,0,0};
    mwSize dimssz[] = {0,0,0};
    mwSize dimssrz[] = {0,0,0};
    mwSize dimss2rz[] = {0,0,0,0};
    mwSize dimsssigmay[] = {0,0,0};
    mwSize dimsssigmaz[] = {0,0,0};


    /* Return rdata structure */
    rdata = (ReturnData*) mxMalloc(sizeof *rdata);
    if (rdata == NULL) return(NULL);

    llhdata = sllhdata = s2llhdata = chi2data = numstepsdata = numrhsevalsdata = orderdata = numstepsSdata =
            numrhsevalsSdata = rzdata = zdata = xdata = ydata = srzdata = szdata = sxdata = sydata = s2rzdata =
            sigmaydata = ssigmaydata = sigmazdata = ssigmazdata = xdotdata = Jdata = dydpdata = dydxdata =
            dxdotdpdata = tsdata = 0;

    mxsol = mxCreateStructMatrix(1,1,29,field_names_sol);

    plhs[0] = mxsol;


    mxstatus = mxCreateDoubleMatrix(1,1,mxREAL);

    mxSetPr(mxstatus,pstatus);
    mxSetField(mxsol,0,"status",mxstatus);

    initField2(llh,1,1);
    initField2(chi2,1,1);
    /*initField2(g,ng,1);
     initField2(r,ng,1);*/

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
        initField2(rz,nmaxevent,nz);
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
                initField3(srz,nmaxevent,nz,np);
                if(sensi>1){
                    initField4(s2rz,nmaxevent,nz,np,np);
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
            if (ng>1) {
                initField2(s2llh,np,(ng-1));
            }
        }
    }

    return(rdata);
}
#endif

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

#ifndef AMICI_WITHOUT_MATLAB
ExpData *setupExpData(const mxArray *prhs[], UserData *udata) {
    /**
     * setupExpData initialises the experimental data struct
     * @param[in] prhs user input @type *mxArray
     * @param[in] udata pointer to the user data struct @type UserData
     * @return edata: experimental data struct @type ExpData
     */

    int nmyt = 0, nmyy = 0, nysigmat = 0, nysigmay = 0; /* integers with problem dimensionality */
    int nmzt = 0, nmzy = 0, nzsigmat = 0, nzsigmay = 0; /* integers with problem dimensionality */

    char *errmsg;
    errmsg = new char[200]();

    ExpData *edata; /* returned rdata struct */


    /* Return edata structure */
    edata = new ExpData();
    if (edata == NULL) return(NULL);

    if ((mxGetM(prhs[8]) == 0 && mxGetN(prhs[8]) == 0) || !prhs[8]) {
        b_expdata = FALSE;
        if(sensi>0 && sensi_meth == AMI_ASA) {
            errMsgIdAndTxt("AMICI:mex:data","No data provied!");
            return NULL;
        }
    } else {
        b_expdata = TRUE;
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

    return(edata);
}
#endif
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void *setupAMI(int *status, UserData *udata, TempData *tdata) {
    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block
     */
    void *ami_mem; /* pointer to ami memory block */
    bool error_corr = TRUE;

    t = tstart;

    g = new realtype[ng]();
    r = new realtype[ng]();

    if (nx > 0) {
        /* allocate temporary objects */
        x = N_VNew_Serial(nx);
        x_old = N_VNew_Serial(nx);
        dx = N_VNew_Serial(nx); /* only needed for idas */
        dx_old = N_VNew_Serial(nx); /* only needed for idas */
        xdot = N_VNew_Serial(nx);
        xdot_old = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);

        if(ne>0) rootsfound = new int[ne]();
        if(ne>0) rootvals= new realtype[ne]();
        if(ne>0) rootidx = new int[nmaxevent*ne*ne]();
        if(ne>0) nroots = new int[ne]();
        if(ne>0) discs = new realtype[nmaxevent*ne]();
        if(ne>0) h = new realtype[ne]();
        if(ne>0) h_tmp = new realtype[ne]();

        if(ne>0) deltax = new realtype[nx]();
        if(ne>0) deltasx = new realtype[nx*np]();
        if(ne>0) deltaxB = new realtype[nx]();
        if(ne>0) deltaqB = new realtype[ng*np]();

        if(ny>0) sigma_y = new realtype[ny]();
        if(ne>0) sigma_z = new realtype[nz]();


        /* initialise states */
        if (x == NULL) return(NULL);
        if(!b_x0) {
            *status = fx0(x, udata);
            if (*status != AMI_SUCCESS) return(NULL);
        } else {
            int ix;
            x_tmp = NV_DATA_S(x);
            for (ix=0; ix<nx; ix++) {
                x_tmp[ix] = x0data[ix];
            }
        }
        *status = fdx0(x, dx, udata); /* only needed for idas */
        if (*status != AMI_SUCCESS) return(NULL);

        /* initialise heaviside variables */
        initHeaviside(status,udata,tdata);

    }

    /* Create AMIS object */
    if (lmm>2||lmm<1) {
        errMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        errMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
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
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* *status = CVLapackDense(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;

             *status = wrap_SetDenseJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;

             break;*/

        case AMI_LAPACKBAND:

            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
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
            *status = AMIKLU(ami_mem, nx, nnz, CSC_MAT);
            if (*status != AMI_SUCCESS) return(NULL);

            *status = wrap_SetSparseJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);

            *status = AMIKLUSetOrdering(ami_mem, ordering);
            if (*status != AMI_SUCCESS) return(NULL);

            break;

        default:
            errMsgIdAndTxt("AMICI:mex:solver","Invalid choice of solver!");
            break;
    }

    if ( sensi >= 1) {

        dydx = new realtype[ny*nx]();
        dydp = new realtype[ny*np]();
        dgdp = new realtype[ng*np*nytrue]();
        dgdx = new realtype[ng*nxtrue*nt]();
        dgdy = new realtype[nytrue*ng*ny]();
        if (ne > 0) {
            dzdp = new realtype[nz*np]();
            dzdx = new realtype[nz*nx]();
        }
        drdp = new realtype[ng*np*nztrue*nmaxevent]();
        drdx = new realtype[ng*nx*nztrue*nmaxevent]();

        dsigma_ydp = new realtype[ny*np]();
        if(ne>0) dsigma_zdp = new realtype[nz*np]();

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
                    int ip;
                    for (ip=0; ip<np; ip++) {
                        sx_tmp = NV_DATA_S(NVsx[plist[ip]]);
                        int ix;
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

                *status = AMIAdjInit(ami_mem, maxsteps, interpType);
                if (*status != AMI_SUCCESS) return(NULL);

                llhS0 = new realtype[ng*np]();
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

void setupAMIB(int *status,void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block for the backward problem
     */
    int ix;

    xB = N_VNew_Serial(nx);
    xB_old = N_VNew_Serial(nx);

    dxB = N_VNew_Serial(nx);

    xQB = N_VNew_Serial(ng*np);
    xQB_old = N_VNew_Serial(ng*np);

    /* write initial conditions */
    if (xB == NULL) return;
    xB_tmp = NV_DATA_S(xB);
    memset(xB_tmp,0,sizeof(realtype)*nx);
    for (ix=0; ix<nx; ix++) {
        xB_tmp[ix] += dgdx[nt-1+ix*nt];
    }
    /*for (ix=0; ix<nxtrue; ix++) {
     for (ig=0; ig<ng; ig++) {
     xB_tmp[ix+ig*nxtrue] += dgdx[nt-1+ix*nt+ig*nxtrue*nt];
     }
     }*/

    if (dxB == NULL) return;
    dxB_tmp = NV_DATA_S(dxB);
    memset(dxB_tmp,0,sizeof(realtype)*nx);

    if (xQB == NULL) return;
    xQB_tmp = NV_DATA_S(xQB);
    memset(xQB_tmp,0,sizeof(realtype)*ng*np);

    /* create backward problem */
    if (lmm>2||lmm<1) {
        errMsgIdAndTxt("AMICI:mex:lmm","Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        errMsgIdAndTxt("AMICI:mex:iter","Illegal value for iter!");
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
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
            break;

        case AMI_LAPACKBAND:


            /* #if SUNDIALS_BLAS_LAPACK
             *status = CVLapackBandB(ami_mem, which, nx, ubw, lbw);
             if (*status != AMI_SUCCESS) return;

             *status = wrap_SetBandJacFnB(ami_mem, which);
             if (*status != AMI_SUCCESS) return;
             #else*/
            errMsgIdAndTxt("AMICI:mex:lapack","Solver currently not supported!");
            /* #endif*/
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
            *status = AMIKLUB(ami_mem, which, nx, nnz, CSC_MAT);
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

void getDataSensisFSA(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getDataSensisFSA extracts data information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ip;
    int iy;
    int ix;

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

    for (iy=0; iy<nytrue; iy++) {
        if(b_expdata){
            if (amiIsNaN(ysigma[iy*nt+it])) {
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
        } else {
            for (ip=0; ip<np; ip++) {
                ssigmaydata[it + nt*(ip*ny+iy)] = amiGetNaN();
            }
        }
    }
    fsy(ts[it],it,sydata,dydx,dydp,NVsx,udata);
    if(b_expdata) {
        fsJy(ts[it],it,sllhdata,s2llhdata,dgdy,dgdp,sydata,dydp,my,udata);
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void prepDataSensis(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * prepDataSensis preprocesses the provided experimental data to compute sensitivities via adjoint or forward methods later on
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int iy,ip,ig;


    *status = fdydx(ts[it],it,dydx,x,udata);
    if (*status != AMI_SUCCESS) return;
    *status = fdydp(ts[it],it,dydp,x,udata);
    if (*status != AMI_SUCCESS) return;
    if(b_expdata) {
        for (iy=0; iy<nytrue; iy++) {
            if (amiIsNaN(ysigma[iy*nt+it])) {
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


        if (sensi_meth == AMI_ASA) {
            for(ig=0; ig<ng; ig++) {
                for(ip=0; ip < np; ip++) {
                    for(iy=0; iy < nytrue; iy++) {
                        if(ig==0) {
                            if (ny>0) {
                                sllhdata[ip] -= dgdp[iy + ip*nytrue];
                            }
                        } else {
                            if (ny>0) {
                                s2llhdata[(ig-1)*np + ip] -= dgdp[(ig*np + ip)*nytrue + iy];
                            }
                        }
                    }
                }
            }
        }
        fdJydy(ts[it],it,dgdy,ydata,x,my,sigma_y,udata);
        fdJydx(ts[it],it,dgdx,ydata,x,dydx,my,sigma_y,udata);
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDataOutput(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getDataOutput extracts output information for data-points
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of current timepoint @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int iy;


    *status = fy(ts[it],it,ydata,x,udata);
    if (*status != AMI_SUCCESS) return;

    if(b_expdata) {
        for (iy=0; iy<nytrue; iy++) {
            /* extract the value for the standard deviation, if the data value is NaN, use
             the parameter value. Store this value in the return struct */
            if (amiIsNaN(ysigma[iy*nt+it])) {
                *status =fsigma_y(t,sigma_y,udata);
                if (*status != AMI_SUCCESS) return;

            } else {
                sigma_y[iy] = ysigma[iy*nt+it];
            }
            sigmaydata[iy*nt+it] = sigma_y[iy];
        }
        fJy(t,it,g,ydata,x,my,sigma_y,udata);
    }
    if (sensi >= 1) {
        prepDataSensis(status, it, ami_mem, udata, rdata, edata, tdata);
        if (sensi_meth == AMI_FSA) {
            getDataSensisFSA(status, it, ami_mem, udata, rdata, edata, tdata);
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisFSA(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getEventSensisFSA extracts event information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */


    *status = fsz(t,ie,nroots,szdata,x,NVsx,udata);
    if (*status != AMI_SUCCESS) return;

}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisFSA_tf(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getEventSensisFSA_tf extracts event information for forward sensitivity
     *     analysis for events that happen at the end of the considered interval
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */


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

void getEventSensisASA(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getEventSensisASA extracts event information for adjoint sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie index of event type @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    int ip;
    int iz;


    for (iz=0; iz<nztrue; iz++) {
        if( z2event[iz]-1 == ie ){
            if(!amiIsNaN(mz[iz*nmaxevent+nroots[ie]])) {
                *status = fdzdp(t,ie,dzdp,x,udata);
                if (*status != AMI_SUCCESS) return;
                *status = fdzdx(t,ie,dzdx,x,udata);
                if (*status != AMI_SUCCESS) return;
                /* extract the value for the standard deviation, if the data value is NaN, use
                 the parameter value. Store this value in the return struct */
                if (amiIsNaN(zsigma[nroots[ie] + nmaxevent*iz])) {
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

                fdJzdp(t,ie,drdp,zdata,x,dzdp,mz,sigma_z,dsigma_zdp,udata,tdata);
                fdJzdx(t,ie,drdx,zdata,x,dzdx,mz,sigma_z,udata,tdata);
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSigma(int *status, int ie, int iz, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getEventSigma extracts fills sigma_z either from the user defined function or from user input
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie event type index @type int
     * @param[in] iz event output index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    /* extract the value for the standard deviation, if the data value is NaN, use
     the parameter value. Store this value in the return struct */
    if (amiIsNaN(zsigma[nroots[ie] + nmaxevent*iz])) {
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

void getEventObjective(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getEventObjective updates the objective function on the occurence of an event
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ie event type index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    if(b_expdata) {
        int iz;
        for (iz=0; iz<nztrue; iz++) {
            if(z2event[iz]-1 == ie) {
                getEventSigma(status, ie, iz, ami_mem, udata, rdata, edata, tdata);
                if(!amiIsNaN(mz[iz*nmaxevent+nroots[ie]])) {
                    r[0] += 0.5*log(2*pi*pow(zsigma[nroots[ie] + nmaxevent*iz],2)) + 0.5*pow( ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] )/zsigma[iz] , 2);
                    *chi2data += pow( ( zdata[nroots[ie] + nmaxevent*iz] - mz[nroots[ie] + nmaxevent*iz] )/zsigma[iz] , 2);
                }
            }
        }
    }

}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventOutput(int *status, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * getEventOutput extracts output information for events
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] tlastroot timepoint of last occured event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    int iz;
    int ie;



    /* EVENT OUTPUT */
    for (ie=0; ie<ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
        if (nroots[ie]<nmaxevent) {
            if(rootsfound[ie] == 1) { /* only consider transitions false -> true */
                *status = fz(t,ie,nroots,zdata,x,udata);
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

                if(b_expdata) {
                    for (iz=0; iz<nztrue; iz++) {
                        if(z2event[iz]-1 == ie) {
                            getEventSigma(status, ie, iz, ami_mem,udata,rdata,edata,tdata);
                            if (*status != AMI_SUCCESS) return;
                        }
                    }

                    getEventObjective(status, ie, ami_mem, udata, rdata, edata, tdata);
                    if (*status != AMI_SUCCESS) return;
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

void fillEventOutput(int *status, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * fillEventOutput fills missing roots at last timepoint
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ie,iz;


    froot(t,x,dx,rootvals,udata);


    /* EVENT OUTPUT */
    if (nztrue>0) {
        for (ie=0; ie<ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
            while (nroots[ie]<nmaxevent) {
                *status = fz(t,ie,nroots,zdata,x,udata);
                if (*status != AMI_SUCCESS) return;


                for (iz=0; iz<nztrue; iz++) {
                    if(z2event[iz]-1 == ie) {
                        rzdata[nroots[ie] + nmaxevent*iz] = rootvals[ie];
                    }
                }


                getEventObjective(status, ie, ami_mem, udata, rdata, edata, tdata);
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

void handleDataPoint(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ix;



    tsdata[it] = ts[it];
    if (nx>0) {
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
    }

    getDataOutput(status, it, ami_mem, udata, rdata, edata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleDataPointB(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * handleDataPoint executes everything necessary for the handling of data points for the backward problems
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it index of data point @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ix;

    xB_tmp = NV_DATA_S(xB);
    for (ix=0; ix<nx; ix++) {
        xB_tmp[ix] += dgdx[it+ix*nt];
    }
    getDiagnosisB(status,it,ami_mem,udata,rdata,tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void handleEvent(int *status, int *iroot, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, ExpData *edata, TempData *tdata, int seflag) {
    /**
     * handleEvent executes everything necessary for the handling of events
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] iroot index of event @type int
     * @param[out] tlastroot pointer to the timepoint of the last event @type *realtype
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[in] seflag flag indicating whether this is a secondary event @type int
     * @return void
     */
    int ie;
    int secondevent = 0;


    /* store heaviside information at event occurence */
    froot(t,x,dx,rootvals,udata);
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
            warnMsgIdAndTxt("AMICI:mex:STUCK_EVENT","AMICI is stuck in an event, as the initial step-size after the event is too small. To fix this, increase absolute and relative tolerances!");
            *status = -99;
            return;
        }
        *tlastroot = t;
    }

    getEventOutput(status, tlastroot, ami_mem, udata, rdata, edata, tdata);
    if (*status != AMI_SUCCESS) return;

    /* if we need to do forward sensitivities later on we need to store the old x and the old xdot */
    if(sensi >= 1){
        /* store x and xdot to compute jump in sensitivities */
        N_VScale(1.0,x,x_old);
        if (sensi_meth == AMI_FSA) {
            *status = fxdot(t,x,dx,xdot,udata);
            N_VScale(1.0,xdot,xdot_old);
            N_VScale(1.0,dx,dx_old);

            /* compute event-time derivative only for primary events, we get into trouble with multiple simultaneously firing events here (but is this really well defined then?), in that case just use the last ie and hope for the best. */
            if (seflag == 0) {
                for (ie = 0; ie<ne; ie++) {
                    if(rootsfound[ie] == 1) { /* only consider transitions false -> true */
                        fstau(t,ie,stau_tmp,x,NVsx,udata);
                    }
                }
            }
        }

        if (sensi_meth == AMI_ASA) {
            /* store x to compute jump in discontinuity */
            if (*iroot<nmaxevent*ne) {
                N_VScale(1.0,x,x_disc[*iroot]);
                N_VScale(1.0,xdot,xdot_disc[*iroot]);
                N_VScale(1.0,xdot_old,xdot_old_disc[*iroot]);
            }
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
        warnMsgIdAndTxt("AMICI:mex:TOO_MUCH_EVENT","Event was recorded but not reported as the number of occured events exceeded (nmaxevents)*(number of events in model definition)!");
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
    froot(t,x,dx,rootvals,udata);
    for (ie = 0; ie<ne; ie++) {
        /* the same event should not trigger itself */
        if (rootsfound[ie] == 0 ) {
            /* check whether there was a zero-crossing */
            if( 0 > h_tmp[ie]*rootvals[ie]) {
                if (h_tmp[ie]<rootvals[ie]) {
                    rootsfound[ie] = 1;
                } else {
                    rootsfound[ie] = -1;
                }
                secondevent++;
            } else {
                rootsfound[ie] = 0;
            }
        } else {
            /* don't fire the same event again */
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

    /* I don't get the meaning of those double if's ... */
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

void handleEventB(int *status, int iroot, void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * handleEventB executes everything necessary for the handling of events for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] iroot index of event @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return cv_status updated status flag @type int
     */

    int ie;
    int ix;
    int ip;
    int ig;


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

            for (ig=0; ig<ng; ig++) {
                for (ip=0; ip<np; ip++) {
                    xQB_tmp[ig*np+ip] += deltaqB[ig*np+ip];
                }
            }


            nroots[ie]--;
        }
    }

    updateHeavisideB(status, iroot, udata, tdata);
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, UserData *udata) {
    /**
     * getTnext computes the next timepoint to integrate to. This is the maximum of
     * tdata and troot but also takes into account if it<0 or iroot<0 where these expressions
     * do not necessarily make sense
     *
     * @param[in] troot timepoint of next event @type realtype
     * @param[in] iroot index of next event @type int
     * @param[in] tdata timepoint of next data point @type realtype
     * @param[in] it index of next data point @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @return tnext next timepoint @type realtype
     */

    realtype tnext;


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

void applyEventBolus(int *status, void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * applyEventBolus applies the event bolus to the current state
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ix;
    int ie;


    for (ie=0; ie<ne; ie++){
        if(rootsfound[ie] == 1) { /* only consider transitions false -> true */
            *status = fdeltax(t,ie,deltax,x,xdot,xdot_old,udata);

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

void applyEventSensiBolusFSA(int *status, void *ami_mem, UserData *udata, TempData *tdata) {
    /**
     * applyEventSensiBolusFSA applies the event bolus to the current sensitivities
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ix;
    int ip;
    int ie;


    for (ie=0; ie<ne; ie++){
        if(rootsfound[ie] == 1) { /* only consider transitions false -> true */
            *status = fdeltasx(t,ie,deltasx,x_old,xdot,xdot_old,NVsx,udata);

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

void initHeaviside(int *status, UserData *udata, TempData *tdata) {
    /**
     * initHeaviside initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ie;


    froot(t,x,dx,rootvals,udata);

    for (ie = 0; ie<ne; ie++) {
        if (rootvals[ie]<=0) {
            h[ie] = 0.0;
        } else {
            h[ie] = 1.0;
        }
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void updateHeaviside(int *status, UserData *udata, TempData *tdata) {
    /**
     * updateHeaviside updates the heaviside variables h on event occurences
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ie;


    /* rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */

    for (ie = 0; ie<ne; ie++) {
        h[ie] += rootsfound[ie];
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void updateHeavisideB(int *status, int iroot, UserData *udata, TempData *tdata) {
    /**
     * updateHeavisideB updates the heaviside variables h on event occurences for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] iroot discontinuity occurance index @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */

    int ie;


    /* rootsfound provides the direction of the zero-crossing, so adding it will give
     the right update to the heaviside variables */

    for (ie = 0; ie<ne; ie++) {
        h[ie] -= rootidx[iroot*ne + ie];
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getDiagnosis(int *status,int it, void *ami_mem, UserData *udata, ReturnData *rdata) {
    /**
     * getDiagnosis extracts diagnosis information from solver memory block and writes them into the return data struct
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return void
     */
    long int numsteps;
    long int numrhsevals;
    int order;


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

void getDiagnosisB(int *status,int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata) {
    /**
     * getDiagnosisB extracts diagnosis information from solver memory block and writes them into the return data struct for the backward problem
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return void
     */
    long int numsteps;
    long int numrhsevals;

    void *ami_memB;


    ami_memB = AMIGetAdjBmem(ami_mem, which);

    *status = AMIGetNumSteps(ami_memB, &numsteps);
    if (*status != AMI_SUCCESS) return;
    numstepsSdata[it] = (realtype)numsteps;

    *status = AMIGetNumRhsEvals(ami_memB, &numrhsevals);
    if (*status != AMI_SUCCESS) return;
    numrhsevalsSdata[it] = (realtype)numrhsevals;

}

#ifdef AMICI_WITHOUT_MATLAB

void initUserDataFields(UserData *udata, ReturnData *rdata) {

    llhdata = sllhdata = s2llhdata = chi2data = numstepsdata = numrhsevalsdata = orderdata = numstepsSdata =
            numrhsevalsSdata = rzdata = zdata = xdata = ydata = srzdata = szdata = sxdata = sydata = s2rzdata =
            sigmaydata = ssigmaydata = sigmazdata = ssigmazdata = xdotdata = Jdata = dydpdata = dydxdata =
            dxdotdpdata = tsdata = 0;

    size_t dimssx[] = {0,0,0};
    size_t dimssy[] = {0,0,0};
    size_t dimssz[] = {0,0,0};
    size_t dimssrz[] = {0,0,0};
    size_t dimss2rz[] = {0,0,0,0};
    size_t dimssigmay[] = {0,0,0};
    size_t dimssigmaz[] = {0,0,0};
    size_t dimsssigmay[] = {0,0,0};
    size_t dimsssigmaz[] = {0,0,0};

    initField2(llh,1,1);
    initField2(chi2,1,1);

    double *mxts = new double[nt]();
    tsdata = mxts;

    initField2(numsteps,nt,1);
    initField2(numrhsevals,nt,1);
    initField2(order,nt,1);
    if(sensi>0){
        initField2(numstepsS,nt,1);
        initField2(numrhsevalsS,nt,1);
    }
    if((nz>0) & (ne>0)){
        initField2(z,nmaxevent,nz);
        initField2(rz,nmaxevent,nz);
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
                initField3(srz,nmaxevent,nz,np);
                if(sensi>1){
                    initField4(s2rz,nmaxevent,nz,np,np);
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
            initField2(s2llh,ng-1,np);
        }
    }
}
#endif


int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, ExpData *edata, int *status, void *ami_mem, int *iroot) {
    /**
     * workForwardProblem solves the forward problem. if forward sensitivities are enabled this will also compute sensitivies
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] iroot pointer to the current root index, the value pointed to will be increased during the forward solve
     * @return int status flag
     */


    /*******************/
    /* FORWARD PROBLEM */
    /*******************/
    int ix, it;
    int ncheck = 0; /* the number of (internal) checkpoints stored so far */
    realtype tlastroot = 0; /* storage for last found root */

    /* loop over timepoints */
    for (it=0; it < nt; it++) {
        if(sensi_meth == AMI_FSA && sensi >= 1) {
            *status = AMISetStopTime(ami_mem, ts[it]);
        }
        if (*status == 0) {
            /* only integrate if no errors occured */
            if(ts[it] > tstart) {
                while (t<ts[it]) {
                    if(sensi_meth == AMI_ASA && sensi >= 1) {
                        if (nx>0) {
                            *status = AMISolveF(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL, &ncheck);
                        } else {
                            t = ts[it];
                        }
                    } else {
                        if (nx>0) {
                            *status = AMISolve(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL);
                        } else {
                            t = ts[it];
                        }
                    }
                    if (nx>0) {
                        x_tmp = NV_DATA_S(x);
                        if (*status == -22) {
                            /* clustering of roots => turn off rootfinding */
                            AMIRootInit(ami_mem, 0, NULL);
                            *status = 0;
                        }
                        /* integration error occured */
                        if (*status<0) {
                            return *status;
                        }
                        if (*status==AMI_ROOT_RETURN) {
                            handleEvent(status, iroot, &tlastroot, ami_mem, udata, rdata, edata, tdata, 0);
                            if (*status != AMI_SUCCESS) return *status;
                        }
                    }
                }
            }

            handleDataPoint(status, it, ami_mem, udata, rdata, edata, tdata);
            if (*status != AMI_SUCCESS) return *status;


        } else {
            for(ix=0; ix < nx; ix++) xdata[ix*nt+it] = amiGetNaN();
        }
    }

    /* fill events */
    if (ne>0) {
        fillEventOutput(status, ami_mem, udata, rdata, edata, tdata);
    }

    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);

    return 0;
}

int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, ExpData *edata, int *status, void *ami_mem, int *iroot, booleantype *setupBdone) {
    /**
     * workBackwardProblem solves the backward problem. if adjoint sensitivities are enabled this will also compute sensitivies
     * workForwardProblem should be called before this is function is called
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] edata pointer to the experimental data struct @type ExpData
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] iroot pointer to the current root index, the value pointed to will be decreased during the forward solve
     * @return int status flag
     */
    int ix, it;
    int ip;

    double tnext;

    if (nx>0) {
        if (sensi >= 1) {
            if(sensi_meth == AMI_ASA) {
                if(*status == 0) {
                    setupAMIB(status, ami_mem, udata, tdata);
                    *setupBdone = true;

                    it = nt-2;
                    (*iroot)--;
                    while (it>=0 || *iroot>=0) {

                        /* check if next timepoint is a discontinuity or a data-point */
                        tnext = getTnext(discs, *iroot, ts, it, udata);

                        if (tnext<t) {
                            *status = AMISolveB(ami_mem, tnext, AMI_NORMAL);
                            if (*status != AMI_SUCCESS) return *status;


                            *status = AMIGetB(ami_mem, which, &t, xB, dxB);
                            if (*status != AMI_SUCCESS) return *status;
                            *status = AMIGetQuadB(ami_mem, which, &t, xQB);
                            if (*status != AMI_SUCCESS) return *status;
                        }

                        /* handle discontinuity */

                        if(ne>0){
                            if(nmaxevent>0){
                                if((*iroot)>=0){
                                    if (tnext == discs[*iroot]) {
                                        handleEventB(status, *iroot, ami_mem, udata, tdata);
                                        (*iroot)--;
                                    }
                                }
                            }
                        }

                        /* handle data-point */

                        if (tnext == ts[it]) {
                            handleDataPointB(status, it, ami_mem, udata, rdata, tdata);
                            it--;
                        }

                        /* reinit states */
                        *status = AMIReInitB(ami_mem, which, t, xB, dxB);
                        if (*status != AMI_SUCCESS) return *status;

                        *status = AMIQuadReInitB(ami_mem, which, xQB);
                        if (*status != AMI_SUCCESS) return *status;

                        *status = AMICalcICB(ami_mem, which, t, xB, dxB);
                        if (*status != AMI_SUCCESS) return *status;
                    }

                    /* we still need to integrate from first datapoint to tstart */
                    if (t>tstart) {
                        if(*status == 0) {
                            if (nx>0) {
                                /* solve for backward problems */
                                *status = AMISolveB(ami_mem, tstart, AMI_NORMAL);
                                if (*status != AMI_SUCCESS) return *status;

                                *status = AMIGetQuadB(ami_mem, which, &t, xQB);
                                if (*status != AMI_SUCCESS) return *status;
                                *status = AMIGetB(ami_mem, which, &t, xB, dxB);
                                if (*status != AMI_SUCCESS) return *status;
                            }
                        }
                    }

                    /* evaluate initial values */
                    NVsx = N_VCloneVectorArray_Serial(np,x);
                    if (NVsx == NULL) return *status;

                    *status = fx0(x,udata);
                    if (*status != AMI_SUCCESS) return *status;
                    *status = fdx0(x,dx,udata);
                    if (*status != AMI_SUCCESS) return *status;
                    *status = fsx0(NVsx, x, dx, udata);
                    if (*status != AMI_SUCCESS) return *status;

                    if(*status == 0) {

                        xB_tmp = NV_DATA_S(xB);

                        int ig;
                        for (ig=0; ig<ng; ig++) {
                            if (ig==0) {
                                for (ip=0; ip<np; ip++) {
                                    llhS0[ig*np + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(NVsx[ip]);
                                    for (ix = 0; ix < nxtrue; ix++) {
                                        llhS0[ip] = llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                                    }
                                }
                            } else {
                                for (ip=0; ip<np; ip++) {
                                    llhS0[ig*np + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(NVsx[ip]);
                                    for (ix = 0; ix < nxtrue; ix++) {
                                        llhS0[ig*np + ip] = llhS0[ig*np + ip] + xB_tmp[ig*nxtrue + ix] * sx_tmp[ix] + xB_tmp[ix] * sx_tmp[ig*nxtrue + ix];
                                    }
                                }
                            }
                        }

                        xQB_tmp = NV_DATA_S(xQB);

                        for(ig=0; ig<ng; ig++) {
                            for(ip=0; ip < np; ip++) {
                                if (ig==0) {
                                    sllhdata[ip] -=  llhS0[ip] + xQB_tmp[ip];
                                    if (nz>0) {
                                        sllhdata[ip] -= drdp[ip];
                                    }
                                } else {
                                    s2llhdata[(ig-1)*np + ip] -= llhS0[ig*np + ip] + xQB_tmp[ig*np + ip];
                                    if (nz>0) {
                                        s2llhdata[(ig-1)*np + ip] -= drdp[ig*np + ip];
                                    }
                                }
                            }
                        }

                    } else {
                        int ig;
                        for(ig=0; ig<ng; ig++) {
                            for(ip=0; ip < np; ip++) {
                                if (ig==0) {
                                    sllhdata[ip] = amiGetNaN();
                                } else {
                                    s2llhdata[(ig-1)*np + ip] = amiGetNaN();
                                }
                            }
                        }
                    }
                } else {
                    int ig;
                    for(ig=0; ig<ng; ig++) {
                        for(ip=0; ip < np; ip++) {
                            if (ig==0) {
                                sllhdata[ip] = amiGetNaN();
                            } else {
                                s2llhdata[(ig-1)*np + ip] = amiGetNaN();
                            }
                        }
                    }
                }
            }
        }
    }

    /* evaluate likelihood */
    if(b_expdata) {
        *llhdata = - g[0] - r[0];
    } else {
        *llhdata = amiGetNaN();
    }

    return 0;
}

void storeJacobianAndDerivativeInReturnData(UserData *udata, TempData *tdata,  ReturnData *rdata) {
    /**
     * evalues the Jacobian and differential equation right hand side, stores it in tdata and
     and copys it to rdata
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return void
     */

    /* store current Jacobian and derivative */
    if(udata) {
        if(tdata) {
            if(nx>0){
                fxdot(t,x,dx,xdot,udata);
                xdot_tmp = NV_DATA_S(xdot);
                memcpy(xdotdata,xdot_tmp,nx*sizeof(realtype));
            }
        }
    }
    if(udata) {
        if(nx>0) {
            fJ(nx,t,0,x,dx,xdot,Jtmp,udata,NULL,NULL,NULL);
            memcpy(Jdata,Jtmp->data,nx*nx*sizeof(realtype));
        }
    }
}

void freeTempDataAmiMem(UserData *udata, TempData *tdata, void *ami_mem, booleantype setupBdone, int status) {
    /**
     * freeTempDataAmiMem frees all allocated memory in udata, tdata and ami_mem
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @param[in] setupBdone flag indicating whether backward problem was initialized @type booleantyp
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[out] status flag indicating success of execution @type *int
     * @return void
     */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(dx);
        N_VDestroy_Serial(xdot);
        N_VDestroy_Serial(x_old);
        N_VDestroy_Serial(dx_old);
        N_VDestroy_Serial(xdot_old);

        delete[] g;
        delete[] r;

        DestroyMat(Jtmp);
        if (ne>0) {
            if(rootsfound) delete[] rootsfound;
            if(rootvals) delete[] rootvals;
            if(rootidx) delete[] rootidx;
            if(sigma_z) delete[] sigma_z;
            if(nroots) delete[] nroots;
            if(discs) delete[] discs;

            if(deltax) delete[] deltax;
            if(deltasx) delete[] deltasx;
            if(deltaxB) delete[] deltaxB;
            if(deltaqB) delete[] deltaqB;
            if(h_tmp) delete[] h_tmp;
        }

        if(ny>0) {
            if(sigma_y)    delete[] sigma_y;
        }
        if (sensi >= 1) {
            if(dydx) delete[] dydx;
            if(dydp) delete[] dydp;
            if(dgdp) delete[] dgdp;
            if(dgdy) delete[] dgdy;
            if(dgdx) delete[] dgdx;
            if(drdp) delete[] drdp;
            if(drdx) delete[] drdx;
            if (ne>0) {
                if(dzdp) delete[] dzdp;
                if(dzdx) delete[] dzdx;
            }
            if(dsigma_ydp) delete[] dsigma_ydp;
            if (ne>0) {
                if(dsigma_zdp) delete[] dsigma_zdp;
            }
            if (sensi_meth == AMI_FSA) {
                N_VDestroyVectorArray_Serial(NVsx,np);
            }
            if (sensi_meth == AMI_ASA) {
                if(NVsx) {
                    N_VDestroyVectorArray_Serial(NVsx,np);
                }
            }

            if (sensi_meth == AMI_FSA) {
                N_VDestroyVectorArray_Serial(sdx, np);
            }
            if (sensi_meth == AMI_ASA) {

                if(llhS0) delete[] llhS0;
                if(setupBdone) N_VDestroy_Serial(dxB);
                if(setupBdone) N_VDestroy_Serial(xB);
                if(setupBdone) N_VDestroy_Serial(xB_old);
                if(setupBdone) N_VDestroy_Serial(xQB);
                if(setupBdone) N_VDestroy_Serial(xQB_old);
            }
        }
        if(ami_mem) N_VDestroy_Serial(id);
        if(ami_mem) AMIFree(&ami_mem);
    }

    if(tdata) delete tdata;
}

#ifdef AMICI_WITHOUT_MATLAB
ReturnData *initReturnData(UserData *udata, int *pstatus) {
    /**
     * initReturnData initialises a ReturnData struct
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] pstatus flag indicating success of execution @type *int
     * @return rdata initialized return data struct @type ReturnData
     */
    ReturnData *rdata; /* returned rdata struct */

    /* Return rdata structure */
    rdata = new ReturnData();
    if (rdata == NULL)
        return(NULL);

    double dblstatus;
    initUserDataFields(udata, rdata);
    *pstatus = (int) dblstatus;

    return(rdata);
}
#endif

#ifdef AMICI_WITHOUT_MATLAB
ReturnData *getSimulationResults(UserData *udata, ExpData *edata, int *pstatus) {
    /**
     * getSimulationResults runs the forward an backwards simulation and returns results in a ReturnData struct
     *
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] edata pointer to the experimental data struct @type ExpData
     * @param[out] pstatus flag indicating success of execution @type *int
     * @return rdata data struct with simulation results @type ReturnData
     */

    double *originalParams = 0;

    if(udata->am_pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * np);
        memcpy(originalParams, p, sizeof(double) * np);

        unscaleParameters(udata);
    }

    int iroot = 0;
    booleantype setupBdone = false;
    *pstatus = 0;
    int problem;
    ReturnData *rdata;
    TempData *tdata = new TempData();
    void *ami_mem = 0; /* pointer to cvodes memory block */
    if (tdata == NULL) goto freturn;


    if (nx>0) {
        ami_mem = setupAMI(pstatus, udata, tdata);
        if (ami_mem == NULL) goto freturn;
    }

    rdata = initReturnData(udata, pstatus);
    if (rdata == NULL) goto freturn;

    if (nx>0) {
        if (edata == NULL) goto freturn;
    }

    *pstatus = 0;

    problem = workForwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot);
    if(problem)
        goto freturn;


    problem = workBackwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot, &setupBdone);
    if(problem)
        goto freturn;

    applyChainRuleFactorToSimulationResults(udata, rdata);

freturn:
    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);
    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, *pstatus);

    if(originalParams) {
        memcpy(p, originalParams, sizeof(double) * np);
        free(originalParams);
    }
    return rdata;
}
#endif


#ifdef AMICI_WITHOUT_MATLAB
void unscaleParameters(UserData *udata) {
    switch(udata->am_pscale) {
    case AMI_SCALING_LOG10:
        for(int i = 0; i < np; ++i) {
            p[i] = pow(10, p[i]);
        }
        break;
    case AMI_SCALING_LN:
        for(int i = 0; i < np; ++i)
            p[i] = exp(p[i]);
        break;
    }
}

void applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata)
{
    if(udata->am_pscale == AMI_SCALING_NONE)
        return;

    // chain-rule factor: multiplier for am_p
    double coefficient;

    switch(udata->am_pscale) {
    case AMI_SCALING_LOG10:
        coefficient = log(10);
        break;
    case AMI_SCALING_LN:
        coefficient = 1.0;
        break;
    }

    if(sensi > 0) {
        if(rdata->am_sllhdata)
            for(int ip = 0; ip < np; ++ip)
                sllhdata[ip] *= p[ip] * coefficient;

        if(rdata->am_sxdata)
            for(int ip = 0; ip < np; ++ip)
                for(int ix = 0; ix < nx; ++ix)
                    for(int it = 0; it < nt; ++it)
                        sxdata[(ip*nx + ix)*nt + it] *= p[ip] * coefficient;

        if(rdata->am_sydata)
            for(int ip = 0; ip < np; ++ip)
                for(int iy = 0; iy < ny; ++iy)
                    for(int it = 0; it < nt; ++it)
                        sydata[(ip*ny + iy)*nt + it] *= p[ip] * coefficient;

        if(rdata->am_szdata)
            for(int ip = 0; ip < np; ++ip)
                for(int iz = 0; iz < nz; ++iz)
                    for(int ie = 0; ie < nmaxevent; ++ie)
                        szdata[(ip*nz + iz)*nmaxevent + ie] *= p[ip] * coefficient;

        if(rdata->am_ssigmaydata)
            for(int ip = 0; ip < np; ++ip)
                for(int iy = 0; iy < ny; ++iy)
                    for(int it = 0; it < nt; ++it)
                        ssigmaydata[(ip*ny + iy)*nt + it] *= p[ip] * coefficient;

        if(rdata->am_ssigmazdata)
            for(int ip = 0; ip < np; ++ip)
                for(int iz = 0; iz < nz; ++iz)
                    for(int ie = 0; ie < nmaxevent; ++ie)
                        ssigmazdata[(ip*nz + iz)*nmaxevent + ie] *= p[ip] * coefficient;
    }

    if(sensi_meth == AMI_SS) {
        if(rdata->am_dxdotdpdata)
            for(int ip = 0; ip < np; ++ip)
                for(int ix = 0; ix < nx; ++ix)
                    dxdotdpdata[ip*nx + ix] *= p[ip] * coefficient;

        if(rdata->am_dydpdata)
            for(int ip = 0; ip < np; ++ip)
                for(int iy = 0; iy < ny; ++iy)
                    dydpdata[ip*nx + iy] *= p[ip] * coefficient;
    }

    if(sensi == 2) {
        if(rdata->am_s2llhdata) {
            double s2coefficient = coefficient * coefficient;
            for(int ip = 0; ip < np; ++ip) {
                for(int iq = 0; iq < np; ++iq) {
                    s2llhdata[ip] *= p[ip] * s2coefficient;
                    if(iq == ip)
                        s2llhdata[ip] += sllhdata[ip];
                }
            }
        }
    }
}

#endif

void processUserData(UserData *udata) {
    /**
     * processUserData initializes fields of the udata struct
     *
     * @param[out] udata pointer to the user data struct @type UserData
     * @return void
     */
    if (nx>0) {
        /* initialise temporary jacobian storage */
        tmp_J = SparseNewMat(nx,nx,nnz,CSC_MAT);
        M_tmp = new realtype[nx*nx]();
        dfdx_tmp = new realtype[nx*nx]();
    }
    if (sensi>0) {
        /* initialise temporary dxdotdp storage */
        tmp_dxdotdp = new realtype[nx*np]();
    }
    if (ne>0) {
        /* initialise temporary stau storage */
        stau_tmp = new realtype[np]();
    }


    w_tmp = new realtype[nw]();
    dwdx_tmp = new realtype[ndwdx]();
    dwdp_tmp = new realtype[ndwdp]();
}
