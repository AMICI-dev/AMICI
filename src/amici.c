/**
 * @file   amici.c
 * @brief  core routines for integration
 */

/** return value indicating successful execution */
#define AMI_SUCCESS               0

UserData setupUserData(const mxArray *prhs[]) {
    /**
     * @brief setupUserData extracts information from the matlab call and returns the corresponding UserData struct
     * @param[in] prhs: pointer to the array of input arguments @type mxArray
     * @return udata: struct containing all provided user data
     */
    
    UserData udata; /* returned udata struct */
    realtype *plistdata; /* input for plist */
    realtype stldetdata; /* input for stldet */
    
    int ip;
    
    /* User udata structure */
    udata = (UserData) mxMalloc(sizeof *udata);
    if (udata == NULL) return(NULL);
    
    /* time */
    
    if (!prhs[1]) {
        mexErrMsgTxt("No time vector provided!");
    }
    ts = mxGetPr(prhs[1]);
    
    nt = (int) mxGetM(prhs[1]) * mxGetN(prhs[1]);
    
    /* parameters */
    
    if (!prhs[2]) {
        mexErrMsgTxt("No parameter vector provided!");
    }
    p = mxGetPr(prhs[2]);
    
    /* constants */
    
    if (!prhs[3]) {
        mexErrMsgTxt("No constant vector provided!");
    }
    k = mxGetPr(prhs[3]);
    
    if (!prhs[4]) {
        mexErrMsgTxt("No options provided!");
    }
    
    if(mxGetField(prhs[4], 0 ,"nx")) { nx = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nx")); } else { mexErrMsgTxt("Parameter nx not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ny")) { ny = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ny")); } else { mexErrMsgTxt("Parameter ny not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"np")) { np = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"np")); } else { mexErrMsgTxt("Parameter np not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ne")) { ne = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ne")); } else { mexErrMsgTxt("Parameter ne not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nnz")) { nnz = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nnz")); } else { mexErrMsgTxt("Parameter nnz not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nmaxevent")) { nmaxevent = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nmaxevent")); } else { mexErrMsgTxt("Parameter nmaxevent not specified as field in options struct!"); }
    
    
    if(mxGetField(prhs[4], 0 ,"tstart")) { tstart = mxGetScalar(mxGetField(prhs[4], 0 ,"tstart")); } else { mexErrMsgTxt("Parameter tstart not specified as field in options struct!"); }
    
    /* plist */
    if (!prhs[5]) {
        mexErrMsgTxt("No parameter list provided!");
    }
    
    if(prhs[5]) {
        plistdata = mxGetPr(prhs[5]);
    }
    
    plist = mxMalloc(np*sizeof(int));
    for (ip=0; ip<np; ip++) {
        plist[ip] = (int)plistdata[ip];
    }
    
    if(mxGetField(prhs[4], 0 ,"atol")) { atol = mxGetScalar(mxGetField(prhs[4], 0 ,"atol")); } else { mexErrMsgTxt("Parameter atol not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"rtol")) { rtol = mxGetScalar(mxGetField(prhs[4], 0 ,"rtol")); } else { mexErrMsgTxt("Parameter rtol not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"maxsteps")) { maxsteps = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"maxsteps")); } else { mexErrMsgTxt("Parameter maxsteps not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"lmm")) { lmm = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"lmm")); } else {  mexErrMsgTxt("Parameter lmm not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"iter")) { iter = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"iter")); } else { mexErrMsgTxt("Parameter iter not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"interpType"))  { interpType = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"interpType")); } else { mexErrMsgTxt("Parameter interpType not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"linsol")) { linsol = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"linsol")); } else { mexErrMsgTxt("Parameter linsol not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"stldet")) { stldetdata = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"stldet")); } else { mexErrMsgTxt("Parameter stldetdata not specified as field in options struct!"); }
    
    if ((int)stldetdata>0.5) {
        stldet = TRUE;
    } else {
        stldet = FALSE;
    }

    if(mxGetField(prhs[4], 0 ,"id")) { idlist = mxGetData(mxGetField(prhs[4], 0, "id")); } else { mexErrMsgTxt("Parameter id not specified as field in options struct!"); }

    if(mxGetField(prhs[4], 0 ,"sensi")) { sensi = (int) mxGetScalar(mxGetField(prhs[4], 0 ,"sensi")); } else { mexErrMsgTxt("Parameter sensi not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ism")) { ism = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ism")); } else { mexErrMsgTxt("Parameter ism not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"sensi_meth")) { sensi_meth = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"sensi_meth")); } else { mexErrMsgTxt("Parameter sensi_meth not specified as field in options struct!"); }
    
    if (sensi > 0) {
        if (sensi_meth != AMI_ASA && sensi_meth != AMI_FSA) {
            mexErrMsgTxt("Invalid sensi_meth specified as field in options struct!");
        }
    }
    
    
    if(mxGetField(prhs[4], 0 ,"ubw")) { ubw = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ubw")); } else { mexErrMsgTxt("Parameter ubw not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"lbw")) { lbw = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"lbw")); } else { mexErrMsgTxt("Parameter lbw not specified as field in options struct!"); }
    
    
    if(mxGetField(prhs[4], 0 ,"sx0")) { sx0data = mxGetPr(mxGetField(prhs[4], 0 ,"sx0")); b_sx0 = TRUE;} else { b_sx0 = FALSE;}
    if (b_sx0) {
        /* check dimensions */
        if(mxGetN(mxGetField(prhs[4], 0 ,"sx0")) != np) { mexErrMsgTxt("Number of rows in sx0 field does not agree with number of model parameters!"); }
        if(mxGetM(mxGetField(prhs[4], 0 ,"sx0")) != nx) { mexErrMsgTxt("Number of columns in sx0 field does not agree with number of model states!"); }
    }
    
    if(mxGetField(prhs[4], 0 ,"data_model")) { data_model = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"data_model")); } else { mexErrMsgTxt("Parameter data_model not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"event_model")) { event_model = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"event_model")); } else { mexErrMsgTxt("Parameter event_model not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"ordering")) { ordering = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ordering")); } else { mexErrMsgTxt("Parameter ordering not specified as field in options struct!"); }

    
    
    /* pbar */
    if (!prhs[6]) {
        mexErrMsgTxt("No parameter scales provided!");
    }
    
    pbar = mxGetPr(prhs[6]);
    
    /* xscale */
    if (!prhs[7]) {
        mexErrMsgTxt("No state scales provided!");
    }
    
    xbar = mxGetPr(prhs[7]);
    
    if (nx>0) {
        /* initialise temporary jacobian storage */
        tmp_J = NewSparseMat(nx,nx,nnz);
    }
    if (sensi>0) {
        /* initialise temporary jacobian storage */
        tmp_dxdotdp = mxMalloc(nx*np*sizeof(realtype));
    }
    
    return(udata);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

ReturnData setupReturnData(const mxArray *prhs[], void *user_data) {
    /**
     * setupReturnData initialises the return data struct
     * @param[in] prhs user input @type *mxArray
     * @param[in] user_data pointer to the user data struct @type UserData
     * @return rdata: return data struct @type ReturnData
     */
    ReturnData rdata; /* returned rdata struct */
    UserData udata; /** user udata */
    
    /* this casting is necessary to ensure availability of accessor macros */
    udata = (UserData) user_data;
    
    /* Return rdata structure */
    rdata = (ReturnData) mxMalloc(sizeof *rdata);
    if (rdata == NULL) return(NULL);
    
    if(mxGetField(prhs[0], 0 ,"t")) { tsdata = mxGetPr(mxGetField(prhs[0], 0 ,"t")); } else { mexErrMsgTxt("t not specified as field in solution struct!"); }
    
    if(mxGetField(prhs[0], 0 ,"x")) { xdata = mxGetPr(mxGetField(prhs[0], 0 ,"x")); } else { mexErrMsgTxt("x not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"y")) { ydata = mxGetPr(mxGetField(prhs[0], 0 ,"y")); } else { mexErrMsgTxt("y not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"xdot")) { xdotdata = mxGetPr(mxGetField(prhs[0], 0 ,"xdot")); } else { mexErrMsgTxt("xdot not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"dxdotdp")) { dxdotdpdata = mxGetPr(mxGetField(prhs[0], 0 ,"dxdotdp")); } else { mexErrMsgTxt("dxdotdp not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"J")) { Jdata = mxGetPr(mxGetField(prhs[0], 0 ,"J")); } else { mexErrMsgTxt("J not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"dydx")) { dydxdata = mxGetPr(mxGetField(prhs[0], 0 ,"dydx")); } else { mexErrMsgTxt("dydx not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"dydp")) { dydpdata = mxGetPr(mxGetField(prhs[0], 0 ,"dydp")); } else { mexErrMsgTxt("dydp not specified as field in solution struct!"); }
    if (ne>0) {
        if(mxGetField(prhs[0], 0 ,"z")) { zdata = mxGetPr(mxGetField(prhs[0], 0 ,"z")); } else { mexErrMsgTxt("z not specified as field in solution struct!"); }
    }
    
    if(mxGetField(prhs[0], 0 ,"numsteps")) { numstepsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numsteps")); } else { mexErrMsgTxt("numsteps not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevals")) { numrhsevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevals")); } else { mexErrMsgTxt("numrhsevals not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"order")) { orderdata = mxGetPr(mxGetField(prhs[0], 0 ,"order")); } else { mexErrMsgTxt("order not specified as field in solution struct!"); }
    
    if(mxGetField(prhs[0], 0 ,"numstepsS")) { numstepsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numstepsS")); } else { mexErrMsgTxt("numstepsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevalsS")) { numrhsevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevalsS")); } else { mexErrMsgTxt("numrhsevalsS not specified as field in solution struct!"); }
    
    if (sensi >= 1) {
        if (sensi_meth == AMI_FSA) {
            if(mxGetField(prhs[0], 0 ,"yS")) { ySdata = mxGetPr(mxGetField(prhs[0], 0 ,"yS")); } else { mexErrMsgTxt("yS not specified as field in solution struct!"); }
            if (ne>0) {
                if(mxGetField(prhs[0], 0 ,"zS")) { zSdata = mxGetPr(mxGetField(prhs[0], 0 ,"zS")); } else { mexErrMsgTxt("zS not specified as field in solution struct!"); }
            }
            if(mxGetField(prhs[0], 0 ,"xS")) { xSdata = mxGetPr(mxGetField(prhs[0], 0 ,"xS")); } else { mexErrMsgTxt("xS not specified as field in solution struct!"); }
        }
        if (sensi_meth == AMI_ASA) {
            if(mxGetField(prhs[0], 0 ,"llhS")) { llhSdata = mxGetPr(mxGetField(prhs[0], 0 ,"llhS")); } else { mexErrMsgTxt("llhS not specified as field in solution struct!"); }
            if (sensi >= 2) {
                if(mxGetField(prhs[0], 0 ,"llhS2")) { llhS2data = mxGetPr(mxGetField(prhs[0], 0 ,"llhS2")); } else { mexErrMsgTxt("llhS not specified as field in solution struct!"); }
            }
        }
    }
    
    if(mxGetField(prhs[0], 0 ,"llh")) { llhdata = mxGetPr(mxGetField(prhs[0], 0 ,"llh")); } else { mexErrMsgTxt("llh not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"chi2")) { chi2data = mxGetPr(mxGetField(prhs[0], 0 ,"chi2")); } else { mexErrMsgTxt("chi2 not specified as field in solution struct!"); }
    
    return(rdata);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

ExpData setupExpData(const mxArray *prhs[], void *user_data) {
    /**
     * setupExpData initialises the experimental data struct
     * @param[in] prhs user input @type *mxArray
     * @param[in] user_data pointer to the user data struct @type UserData
     * @return edata: experimental data struct @type ExpData
     */
    
    int nmyt, nmyy, nysigmat, nysigmay; /* integers with problem dimensionality */
    int nmzt, nmzy, nzsigmat, nzsigmay; /* integers with problem dimensionality */
    
    ExpData edata; /* returned rdata struct */
    UserData udata; /** user udata */
    
    udata = (UserData) user_data;
    
    /* Return rdata structure */
    edata = (ExpData) mxMalloc(sizeof *edata);
    if (edata == NULL) return(NULL);
    
    if (data_model == AMI_ONEOUTPUT) {
        if ( (ny>1) | (nt>1) ) {
            mexErrMsgTxt("Data model AMI_ONEOUTPUT not allowed for more than one time-point or more than one observable!");
        }
    } else {
        
        if (!prhs[8]) {
            mexErrMsgTxt("No data provided!");
        }
        if (mxGetField(prhs[8], 0 ,"Y")) {
            my = mxGetPr(mxGetField(prhs[8], 0 ,"Y"));
            nmyy = (int) mxGetN(mxGetField(prhs[8], 0 ,"Y"));
            nmyt = (int) mxGetM(mxGetField(prhs[8], 0 ,"Y"));
        } else {
            mexErrMsgTxt("Field Y not specified as field in data struct!");
        }
        
        if (mxGetField(prhs[8], 0 ,"Sigma_Y")) {
            ysigma = mxGetPr(mxGetField(prhs[8], 0 ,"Sigma_Y"));
            nysigmay = (int) mxGetN(mxGetField(prhs[8], 0 ,"Sigma_Y"));
            nysigmat = (int) mxGetM(mxGetField(prhs[8], 0 ,"Sigma_Y"));
        } else {
            mexErrMsgTxt("Field Sigma_Y not specified as field in data struct!");
        }
        if (mxGetField(prhs[8], 0 ,"T")) {
            mz = mxGetPr(mxGetField(prhs[8], 0 ,"T"));
            nmzy = (int) mxGetN(mxGetField(prhs[8], 0 ,"T"));
            nmzt = (int) mxGetM(mxGetField(prhs[8], 0 ,"T"));
        } else {
            mexErrMsgTxt("Field T not specified as field in data struct!");
        }
        
        if (mxGetField(prhs[8], 0 ,"sigma_z")) {
            zsigma = mxGetPr(mxGetField(prhs[8], 0 ,"sigma_z"));
            nzsigmay = (int) mxGetN(mxGetField(prhs[8], 0 ,"sigma_z"));
            nzsigmat = (int) mxGetM(mxGetField(prhs[8], 0 ,"sigma_z"));
        } else {
            mexErrMsgTxt("Field sigma_z not specified as field in data struct!");
        }
        
        if (nmyt != nt) {
            mexErrMsgTxt("Number of time-points in data matrix does not match provided time vector");
        }
        
        if (nysigmat != nt) {
            mexErrMsgTxt("Number of time-points in data-sigma matrix does not match provided time vector");
        }
        
        if (nmyy != ny) {
            mexErrMsgTxt("Number of observables in data matrix does not match provided ny");
        }
        
        if (nysigmay != ny) {
            mexErrMsgTxt("Number of observables in data-sigma matrix does not match provided ny");
        }
        
        if (nmzt != nmaxevent) {
            mexErrMsgTxt("Number of time-points in event matrix does not match provided nmaxevent");
        }
        
        if (nzsigmat != nmaxevent) {
            mexErrMsgTxt("Number of time-points in event-sigma matrix does not match provided nmaxevent");
        }
        
        if (nmzy != ne) {
            mexErrMsgTxt("Number of events in event matrix does not match provided ne");
        }
        
        if (nzsigmay != ne) {
            mexErrMsgTxt("Number of events in event-sigma matrix does not match provided ne");
        }
        
    }
    
    
    return(edata);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void *setupAMI(int *status, void *user_data, void *temp_data) {
    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[in] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    void *ami_mem; /* pointer to ami memory block */
    bool error_corr = TRUE;
    int pretype= PREC_NONE; /* specifies the type of preconditioning and must be one of: PREC NONE, PREC LEFT, PREC RIGHT, or PREC BOTH. */
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    TempData tdata; /* user udata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    int ip;
    int ix;
    
    r = 0;
    g = 0;
    
    if (nx > 0) {
        
        /* write initial conditions */
        x = N_VNew_Serial(nx);
        dx = N_VNew_Serial(nx); /* only needed for idas */
        xdot = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);
        
        if(ne>0) rootvaltmp = mxMalloc(ne*sizeof(realtype));
        if(ne>0) rootsfound = mxMalloc(ne*sizeof(int));
        if(ne>0) rootidx = mxMalloc(nmaxevent*sizeof(int));
        if(ne>0) discs = mxMalloc(nmaxevent*sizeof(realtype));
        
        if(ny>0) sigma_y = mxMalloc(ny*sizeof(realtype));
        if(ny>0) memset(sigma_y,0,ny*sizeof(realtype));
        if(ne>0) sigma_z = mxMalloc(ne*sizeof(realtype));
        if(ne>0) memset(sigma_z,0,ne*sizeof(realtype));
        
        if (x == NULL) return(NULL);
        *status = fx0(x, udata);
        if (*status != AMI_SUCCESS) return(NULL);
        *status = fdx0(x, dx, udata); /* only needed for idas */
        if (*status != AMI_SUCCESS) return(NULL);
    }
    
    /* Create AMIS object */
    if (lmm>2||lmm<1) {
        mexErrMsgTxt("Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        mexErrMsgTxt("Illegal value for iter!");
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
    
    *status = AMISetUserData(ami_mem, udata); /* attaches userdata*/
    if(*status != AMI_SUCCESS) return(NULL);
    
    *status = AMISetMaxNumSteps(ami_mem, maxsteps); /* specify maximal number of steps */
    if(*status != AMI_SUCCESS) return(NULL);
    
    *status = AMISetStabLimDet(ami_mem, stldet); /* activates stability limit detection */
    if(*status != AMI_SUCCESS) return(NULL);
    
    if (ne>0) {
        *status = wrap_RootInit(ami_mem, udata); /* activates root detection */
        if(*status != AMI_SUCCESS) return(NULL);
    }
    
    /* Attach linear solver module */
    switch (linsol) {
            
            /* DieECT SOLVERS */
            
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
            mexErrMsgTxt("Solver currently not supported!");
            /* *status = CVLapackDense(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetDenseJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;
             
             break;*/
            
        case AMI_LAPACKBAND:
            
            mexErrMsgTxt("Solver currently not supported!");
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
            
        case AMI_KLU:
            *status = AMIKLU(ami_mem, nx, nnz);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetSparseJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
        
            *status = AMIKLUSetOrdering(ami_mem, ordering);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        default:
            mexErrMsgTxt("Invalid choice of solver!");
            break;
    }
    
    if ( sensi >= 1) {
        
        if (sensi_meth == AMI_FSA) {
            
            if(nx>0) {
                
                /* Set sensitivity initial conditions */
                
                sx = N_VCloneVectorArray_Serial(np, x);
                sdx = N_VCloneVectorArray_Serial(np, x);
                if (sx == NULL) return(NULL);
                if (sdx == NULL) return(NULL);
                
                if(!b_sx0) {
                    *status = fsx0(sx, x, dx, udata);
                    if (*status != AMI_SUCCESS) return(NULL);
                } else {
                    for (ip=0; ip<np; ip++) {
                        sx_tmp = NV_DATA_S(sx[plist[ip]]);
                        for (ix=0; ix<nx; ix++) {
                            sx_tmp[ix] = sx0data[ix + nx*plist[ip]];
                        }
                    }
                }
                *status = fsdx0(sdx, x, dx, udata);
                if (*status != AMI_SUCCESS) return(NULL);
                
                /* Activate sensitivity calculations */
                
                *status = wrap_SensInit1(ami_mem, sx, sdx, udata);
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
                
                *status = AMIAdjInit(ami_mem, maxsteps, interpType);
                if (*status != AMI_SUCCESS) return(NULL);
                
                dydx = mxMalloc(ny*nx*sizeof(realtype));
                memset(dydx,0,ny*nx*sizeof(realtype));
                dydp = mxMalloc(ny*np*sizeof(realtype));
                memset(dydp,0,ny*np*sizeof(realtype));
                llhS0 = mxMalloc(np*sizeof(realtype));
                memset(llhS0,0,np*sizeof(realtype));
                dgdp = mxMalloc(np*sizeof(realtype));
                memset(dgdp,0,np*sizeof(realtype));
                dgdx = mxMalloc(nx*nt*sizeof(realtype));
                memset(dgdx,0,nx*nt*sizeof(realtype));
                if (ne > 0) {
                    dzdp = mxMalloc(ne*np*sizeof(realtype));
                    memset(dzdp,0,ne*np*sizeof(realtype));
                    dzdx = mxMalloc(ne*nx*sizeof(realtype));
                    memset(dzdx,0,ne*nx*sizeof(realtype));
                }
                drdp = mxMalloc(np*sizeof(realtype));
                memset(drdp,0,np*sizeof(realtype));
                drdx = mxMalloc(nx*nmaxevent*sizeof(realtype));
                memset(drdx,0,nx*nmaxevent*sizeof(realtype));
                dsigma_ydp = mxMalloc(ny*np*sizeof(realtype));
                memset(dsigma_ydp,0,ny*np*sizeof(realtype));
                if(ne>0) dsigma_zdp = mxMalloc(ne*np*sizeof(realtype));
                if(ne>0) memset(dsigma_zdp,0,ne*np*sizeof(realtype));
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

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void setupAMIB(int *status,void *ami_mem, void *user_data, void *temp_data) {
    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type *int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[in] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    bool error_corr = TRUE;
    int pretype= PREC_NONE; /* specifies the type of preconditioning and must be one of: PREC NONE, PREC LEFT, PREC RIGHT, or PREC BOTH. */
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    tdata = (TempData) temp_data;
    
    int ix;
    
    xB = N_VNew_Serial(nx);
    
    dxB = N_VNew_Serial(nx);
    
    xQB = N_VNew_Serial(np);
    
    /* BACKWARD PROBLEM */
    
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
        mexErrMsgTxt("Illegal value for lmm!");
    }
    if (iter>2||iter<1) {
        mexErrMsgTxt("Illegal value for iter!");
    }
    *status = AMICreateB(ami_mem, lmm, iter, &which);
    if (*status != AMI_SUCCESS) return;
    
    
    /* allocate memory for the backward problem */
    *status = wrap_binit(ami_mem, which, xB, dxB, t);
    if (*status != AMI_SUCCESS) return;
    
    /* specify integration tolerances for backward problem */
    *status = AMISStolerancesB(ami_mem, which, RCONST(rtol), RCONST(atol));
    if(*status != AMI_SUCCESS) return;
    
    /* Attach user data */
    *status = AMISetUserDataB(ami_mem, which, udata);
    if(*status != AMI_SUCCESS) return;
    
    /* Number of maximal internal steps */
    *status = AMISetMaxNumStepsB(ami_mem, which, maxsteps);
    if(*status != AMI_SUCCESS) return;
    
    switch (linsol) {
            
            /* DieECT SOLVERS */
            
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
            mexErrMsgTxt("Solver currently not supported!");
            /* #endif*/
            break;
            
        case AMI_LAPACKBAND:
            
            
            /* #if SUNDIALS_BLAS_LAPACK
             *status = CVLapackBandB(ami_mem, which, nx, ubw, lbw);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetBandJacFnB(ami_mem, which);
             if (*status != AMI_SUCCESS) return;
             #else*/
            mexErrMsgTxt("Solver currently not supported!");
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

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

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
        
        if(data_model != AMI_ONEOUTPUT) {
            if (mxIsNaN(ysigma[iy*nt+it])) {
                *status =fsigma_y(t,sigma_y,udata);
                if (*status != AMI_SUCCESS) return;
                
            } else {
                sigma_y[iy] = ysigma[iy*nt+it];
            }
        }
        
        if (data_model == AMI_NORMAL) {
            if(!mxIsNaN(my[iy*nt+it])){
                g += 0.5*log(2*pi*pow(sigma_y[iy],2)) + 0.5*pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/sigma_y[iy] , 2);
                *chi2data += pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/sigma_y[iy] , 2);
            }
        }
        if (data_model == AMI_LOGNORMAL) {
            if(!mxIsNaN(my[iy*nt+it])){
                g += 0.5*log(2*pi*pow(sigma_y[iy]*ydata[iy*nt+it],2)) + 0.5*pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/sigma_y[iy] , 2);
                *chi2data += pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/sigma_y[iy] , 2);
            }
        }
        if (data_model == AMI_ONEOUTPUT) {
            g += ydata[iy*nt+it];
        }
    }
    if (sensi >= 1) {
        if (sensi_meth == AMI_FSA) {
            getDataSensisFSA(&status, it, ami_mem, udata, rdata, edata, tdata);
        }
        if (sensi_meth == AMI_ASA) {
            getDataSensisASA(&status, it, ami_mem, udata, rdata, edata, tdata);
        }
    }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int getEventOutput(int *status, realtype *tlastroot, int *nroots, int *iroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventOutput extracts output information for events
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
    
    int cv_status;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    if (t == *tlastroot) {
        /* we are stuck in a root => turn off rootfinding */
        /* at some point we should find a more intelligent solution here, and turn on rootfinding again after some time */
        AMIRootInit(ami_mem, 0, NULL);
        cv_status = 0;
    }
    *tlastroot = t;
    if (sensi >= 1) {
        if(sensi_meth == AMI_ASA) {
            getEventSensisASA(status, nroots, iroot, ami_mem, udata, rdata, edata, tdata);
        } else {
            getEventSensisFSA(status, nroots, ami_mem, udata, rdata, tdata);
        }
    }
    if (cv_status == -22) {
        /* clustering of roots => turn off rootfinding */
        AMIRootInit(ami_mem, 0, NULL);
        cv_status = 0;
    }
    return(cv_status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void fillEventOutput(int *status, int *nroots, int *iroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventOutput fills missing roots with last timepoint
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
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;

    if (sensi >= 1) {
        if(sensi_meth == AMI_ASA) {
            getEventSensisASA(status, nroots, iroot, ami_mem, udata, rdata, edata, tdata);
        } else {
            getEventSensisFSA(status, nroots, ami_mem, udata, rdata, tdata);
        }
    }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

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
                *status = AMIGetSens(ami_mem, &t, sx);
                if (*status != AMI_SUCCESS) return;
            }
            
            sx_tmp = NV_DATA_S(sx[ip]);
            for(ix=0; ix < nx; ix++) {
                xSdata[(ip*nx + ix)*nt + it] = sx_tmp[ix];
            }
        }
    }
    fsy(ts[it],it,ySdata,x,sx,udata);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

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
    int ix;
    
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
        if(data_model != AMI_ONEOUTPUT) {
            if (mxIsNaN(ysigma[iy*nt+it])) {
                *status = fsigma_y(t,sigma_y,udata);
                if (*status != AMI_SUCCESS) return;
                *status = fdsigma_ydp(t,dsigma_ydp,udata);
                if (*status != AMI_SUCCESS) return;
            } else {
                for (ip=0; ip<np; ip++) {
                    dsigma_ydp[ip*ny+iy] = 0;
                }
                sigma_y[iy] = ysigma[iy*nt+it];
            }
        }
        for (ip=0; ip<np; ip++) {
            if(data_model == AMI_NORMAL) {
                if(!mxIsNaN(my[iy*nt+it])){
                    dgdp[ip] += dsigma_ydp[ip*ny+iy]/sigma_y[iy] + ( dydp[ip*ny+iy]* ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( sigma_y[iy] , 2) - dsigma_ydp[ip*ny+iy]*pow( ( ydata[iy*nt+it] - my[iy*nt+it] ),2)/pow( sigma_y[iy] , 3);
                }
            }
            if(data_model == AMI_LOGNORMAL) {
                if(!mxIsNaN(my[iy*nt+it])){
                    dgdp[ip] += (sigma_y[iy]*dydp[ip*ny+iy] + ydata[iy*nt+it]*dsigma_ydp[ip*ny+iy])/(sigma_y[iy]*ydata[iy*nt+it]) + ( dydp[ip*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( sigma_y[iy] , 2) -  dsigma_ydp[ip*ny+iy]*pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ),2)/pow(sigma_y[iy] , 3);
                }
            }
            if(data_model == AMI_ONEOUTPUT) {
                dgdp[ip] += dydp[ip*ny+iy];
            }
        }
        for (ix=0; ix<nx; ix++) {
            if(data_model == AMI_NORMAL) {
                if(!mxIsNaN(my[iy*nt+it])){
                    dgdx[it+ix*nt] += ( dydx[ix*ny+iy] * ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( sigma_y[iy] , 2);
                }
            }
            if(data_model == AMI_LOGNORMAL) {
                if(!mxIsNaN(my[iy*nt+it])){
                    dgdx[it+ix*nt] += 1/(2*pi)*dydx[ix*ny+iy]/ydata[iy*nt+it] + ( dydx[ix*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( sigma_y[iy] , 2);
                }
            }
            if(data_model == AMI_ONEOUTPUT) {
                dgdx[it+ix*nt] += dydx[ix*ny+iy];
            }
        }
    }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void getEventSensisFSA(int *status, int *nroots, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getEventSensisFSA extracts root information for forward sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type int
     * @param[out] nroots counter for the number of found roots @type int
     * @param[in] ami_mem pointer to the solver memory block @type void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    
    int ip;
    int jp;
    int ie;
    
    /* this casting is necessary to ensure availability of accessor macros */
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    TempData tdata; /* temp data */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    tdata = (TempData) temp_data;
    
    *status = AMIGetRootInfo(ami_mem, rootsfound);
    if (*status != AMI_SUCCESS) return;
    /* ROOTS FOR ROOTFUNCTION */
    if (*nroots<nmaxevent) {
        for (ie=0; ie<ne; ie++){ /* only look for roots of the rootfunction not discontinuities */
            if(rootsfound[ie] != 0) {
                *status = fz(t,*nroots,zdata,x,udata);
                /* extract sensitivity information */
                if(sensi >= 1) {
                    if(sensi_meth == AMI_FSA) {
                        *status = AMIGetSens(ami_mem, &t, sx);
                        if (*status != AMI_SUCCESS) return;
                        *status = fsz(t,*nroots,zSdata,x,sx,udata);
                        if (*status != AMI_SUCCESS) return;
                    }
                }
                (*nroots)++;
            }
        }
    }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void getEventSensisASA(int *status, int *nroots, int *iroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getEventSensisASA extracts root information for adjoint sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] nroots counter for the number of found roots @type *int
     * @param[out] iroot counter for the number of found discontinuities @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int ip;
    int ie;
    int ix;
    
    UserData udata; /* user udata */
    ReturnData rdata; /* return rdata */
    ExpData edata; /* exp edata */
    TempData tdata; /* temp tdata */
    udata = (UserData) user_data;
    rdata = (ReturnData) return_data;
    edata = (ExpData) exp_data;
    tdata = (TempData) temp_data;
    
    
    *status = AMIGetRootInfo(ami_mem, rootsfound);
    if (*status != AMI_SUCCESS) return;
    /* EVENT ROOTS */
    if (*nroots<nmaxevent) {
        for (ie=0; ie<ne; ie++){
            if(rootsfound[ie] == 1) {
                *status = fz(t,*nroots,zdata,x,udata);
                /* extract sensitivity information */
                if(!mxIsNaN(mz[ie*nmaxevent+*nroots])) {
                    if (event_model == AMI_NORMAL) {
                        r += 0.5*log(2*pi*pow(zsigma[ie*nmaxevent+*nroots],2)) + 0.5*pow( ( zdata[ie] - mz[ie*nmaxevent+*nroots] )/zsigma[ie*nmaxevent+*nroots] , 2);
                        *chi2data += pow( ( zdata[ie] - mz[ie*nmaxevent+*nroots] )/zsigma[ie*nmaxevent+*nroots] , 2);
                    }
                    if (sensi>=1) {
                        *status = fdzdp(t,ie,dzdp,x,udata);
                        if (*status != AMI_SUCCESS) return;
                        *status = fdzdx(t,ie,dzdx,x,udata);
                        if (*status != AMI_SUCCESS) return;
                        if (mxIsNaN(zsigma[ie*nmaxevent + *nroots])) {
                            *status = fsigma_z(t,ie,sigma_z,udata);
                            if (*status != AMI_SUCCESS) return;
                            *status = fdsigma_zdp(t,ie,dsigma_zdp,udata);
                            if (*status != AMI_SUCCESS) return;
                        } else {
                            for (ip=0; ip<np; ip++) {
                                dsigma_zdp[ip*ne+ie] = 0;
                            }
                            sigma_z[ie] = zsigma[ie*nmaxevent + *nroots];
                        }
                        
                        for (ip=0; ip<np; ip++) {
                            if(event_model == AMI_NORMAL) {
                                drdp[ip] += dsigma_zdp[ip*ne+ie]/sigma_z[ie] + ( dzdp[ip*ne+ie]* ( zdata[ie] - mz[ie*nmaxevent+*nroots] ) )/pow( sigma_z[ie] , 2) - dsigma_zdp[ip*ne+ie]*pow( zdata[ie] - mz[ie*nmaxevent+*nroots] ,2)/pow( sigma_z[ie] , 3);
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(event_model  == AMI_NORMAL) {
                                drdx[*nroots+ix*nmaxevent] += ( dzdx[ix*ne+ie] * ( zdata[ie] - mz[ie*nmaxevent+*nroots] ) )/pow( sigma_z[ie] , 2);
                            }
                        }
                    }
                }
                (*nroots)++;
            }
        }
    }
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

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

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void getDiagnosisB(int *status,int it, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getDiagnosis extracts diagnosis information from solver memory block and writes them into the return data struct
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