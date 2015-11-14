/**
 * @file   amici.c
 * @brief  core routines for integration
 */

#ifndef amici_c
/** include guard */
#define amici_c
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
    if(mxGetField(prhs[4], 0 ,"nr")) { nr = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nr")); } else { mexErrMsgTxt("Parameter nr not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ndisc")) { ndisc = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ndisc")); } else { mexErrMsgTxt("Parameter ndisc not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nnz")) { nnz = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nnz")); } else { mexErrMsgTxt("Parameter nnz not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nmaxroot")) { nmaxroot = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nmaxroot")); } else { mexErrMsgTxt("Parameter nmaxroot not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nmaxdisc")) { nmaxdisc = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nmaxdisc")); } else { mexErrMsgTxt("Parameter nmaxdisc not specified as field in options struct!"); }
    
    
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
        if (sensi_meth != CW_ASA && sensi_meth != CW_FSA) {
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
    rval = 0;
    
    if (nx > 0) {
        
        /* write initial conditions */
        x = N_VNew_Serial(nx);
        dx = N_VNew_Serial(nx); /* only needed for idas */
        xdot = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);
        
        if(nr>0) rootvaltmp = mxMalloc(nr*sizeof(realtype));
        if(nr+ndisc>0) rootsfound = mxMalloc((nr+ndisc)*sizeof(int));
        if(nr>0) rootidx = mxMalloc(nr*sizeof(int));
        if(ndisc>0) discs = mxMalloc(nmaxdisc*sizeof(realtype));
        if(ndisc>0) irdiscs = mxMalloc(nmaxdisc*sizeof(realtype));
        
        if(ny>0) sigma_y = mxMalloc(ny*sizeof(realtype));
        if(ny>0) memset(sigma_y,0,ny*sizeof(realtype));
        if(nr>0) sigma_t = mxMalloc(nr*sizeof(realtype));
        if(nr>0) memset(sigma_t,0,nr*sizeof(realtype));
        
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
    
    if (nr+ndisc>0) {
        *status = wrap_RootInit(ami_mem, udata); /* activates root detection */
        if(*status != AMI_SUCCESS) return(NULL);
    }
    
    /* Attach linear solver module */
    switch (linsol) {
            
            /* DIRECT SOLVERS */
            
        case CW_DENSE:
            *status = AMIDense(ami_mem, nx);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetDenseJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case CW_BAND:
            *status = AMIBand(ami_mem, nx, ubw, lbw);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetBandJacFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case CW_LAPACKDENSE:
            mexErrMsgTxt("Solver currently not supported!");
            /* *status = CVLapackDense(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetDenseJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;
             
             break;*/
            
        case CW_LAPACKBAND:
            
            mexErrMsgTxt("Solver currently not supported!");
            /* *status = CVLapackBand(ami_mem, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetBandJacFn(ami_mem);
             if (*status != AMI_SUCCESS) return;
             
             break;*/
            
        case CW_DIAG:
            *status = AMIDiag(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case CW_SPGMR:
            *status = AMISpgmr(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);

            break;
            
        case CW_SPBCG:
            *status = AMISpbcg(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case CW_SPTFQMR:
            *status = AMISptfqmr(ami_mem, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return(NULL);
            
            *status = wrap_SetJacTimesVecFn(ami_mem);
            if (*status != AMI_SUCCESS) return(NULL);
            
            break;
            
        case CW_KLU:
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
        
        if (sensi_meth == CW_FSA) {
            
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
        
        if (sensi_meth == CW_ASA) {
            
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
                if (nr > 0) {
                    dtdp = mxMalloc(nr*np*sizeof(realtype));
                    memset(dgdp,0,nr*np*sizeof(realtype));
                    dtdx = mxMalloc(nr*nx*sizeof(realtype));
                    memset(dgdx,0,nr*nx*sizeof(realtype));
                }
                drdp = mxMalloc(np*sizeof(realtype));
                memset(drdp,0,np*sizeof(realtype));
                drdx = mxMalloc(nx*nmaxroot*sizeof(realtype));
                memset(drdx,0,nx*nmaxroot*sizeof(realtype));
                drvaldp = mxMalloc(np*sizeof(realtype));
                memset(drvaldp,0,np*sizeof(realtype));
                drvaldx = mxMalloc(nx*nmaxroot*sizeof(realtype));
                memset(drvaldx,0,nx*nmaxroot*sizeof(realtype));
                dsigma_ydp = mxMalloc(ny*np*sizeof(realtype));
                memset(dsigma_ydp,0,ny*np*sizeof(realtype));
                if(nr>0) dsigma_tdp = mxMalloc(nr*np*sizeof(realtype));
                if(nr>0) memset(dsigma_tdp,0,nr*np*sizeof(realtype));
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
            
            /* DIRECT SOLVERS */
            
        case CW_DENSE:
            *status = AMIDenseB(ami_mem, which, nx);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetDenseJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case CW_BAND:
            *status = AMIBandB(ami_mem, which, nx, ubw, lbw);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetBandJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case CW_LAPACKDENSE:
            
            /* #if SUNDIALS_BLAS_LAPACK
             *status = CVLapackDenseB(ami_mem, which, nx);
             if (*status != AMI_SUCCESS) return;
             
             *status = wrap_SetDenseJacFnB(ami_mem, which);
             if (*status != AMI_SUCCESS) return;
             #else*/
            mexErrMsgTxt("Solver currently not supported!");
            /* #endif*/
            break;
            
        case CW_LAPACKBAND:
            
            
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
            
        case CW_DIAG:
            *status = AMIDiagB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetDenseJacFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
            /* ITERATIVE SOLVERS */
            
        case CW_SPGMR:
            *status = AMISpgmrB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case CW_SPBCG:
            *status = AMISpbcgB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case CW_SPTFQMR:
            *status = AMISptfqmrB(ami_mem, which, PREC_NONE, 5);
            if (*status != AMI_SUCCESS) return;
            
            *status = wrap_SetJacTimesVecFnB(ami_mem, which);
            if (*status != AMI_SUCCESS) return;
            
            break;
            
        case CW_KLU:
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
    if (nr>0) {
        if(mxGetField(prhs[0], 0 ,"root")) { rootdata = mxGetPr(mxGetField(prhs[0], 0 ,"root")); } else { mexErrMsgTxt("root not specified as field in solution struct!"); }
        if(mxGetField(prhs[0], 0 ,"rootval")) { rootvaldata = mxGetPr(mxGetField(prhs[0], 0 ,"rootval")); } else { mexErrMsgTxt("root not specified as field in solution struct!"); }
    }
    
    if(mxGetField(prhs[0], 0 ,"numsteps")) { numstepsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numsteps")); } else { mexErrMsgTxt("numsteps not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevals")) { numrhsevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevals")); } else { mexErrMsgTxt("numrhsevals not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"order")) { orderdata = mxGetPr(mxGetField(prhs[0], 0 ,"order")); } else { mexErrMsgTxt("order not specified as field in solution struct!"); }
    
    if(mxGetField(prhs[0], 0 ,"numstepsS")) { numstepsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numstepsS")); } else { mexErrMsgTxt("numstepsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevalsS")) { numrhsevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevalsS")); } else { mexErrMsgTxt("numrhsevalsS not specified as field in solution struct!"); }
    
    if (sensi >= 1) {
        if (sensi_meth == CW_FSA) {
            if(mxGetField(prhs[0], 0 ,"yS")) { ySdata = mxGetPr(mxGetField(prhs[0], 0 ,"yS")); } else { mexErrMsgTxt("yS not specified as field in solution struct!"); }
            if (nr>0) {
                if(mxGetField(prhs[0], 0 ,"rootS")) { rootSdata = mxGetPr(mxGetField(prhs[0], 0 ,"rootS")); } else { mexErrMsgTxt("rootS not specified as field in solution struct!"); }
                if(mxGetField(prhs[0], 0 ,"rootvalS")) { rootvalSdata = mxGetPr(mxGetField(prhs[0], 0 ,"rootvalS")); } else { mexErrMsgTxt("rootvalS not specified as field in solution struct!"); }
                if (sensi >= 2) {
                    if(mxGetField(prhs[0], 0 ,"rootS2")) { rootS2data = mxGetPr(mxGetField(prhs[0], 0 ,"rootS2")); } else { mexErrMsgTxt("rootS2 not specified as field in solution struct!"); }
                    if(mxGetField(prhs[0], 0 ,"rootvalS2")) { rootvalS2data = mxGetPr(mxGetField(prhs[0], 0 ,"rootvalS2")); } else { mexErrMsgTxt("rootvalS2 not specified as field in solution struct!"); }
                }
            }
            if(mxGetField(prhs[0], 0 ,"xS")) { xSdata = mxGetPr(mxGetField(prhs[0], 0 ,"xS")); } else { mexErrMsgTxt("xS not specified as field in solution struct!"); }
        }
        if (sensi_meth == CW_ASA) {
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

ExpData setupExpData(const mxArray *prhs[], void *user_data) {
        /**
     * setupExpData initialises the experimental data struct
     * @param[in] prhs user input @type *mxArray
     * @param[in] user_data pointer to the user data struct @type UserData
     * @return edata: experimental data struct @type ExpData
     */
    
    int nmyt, nmyy, nysigmat, nysigmay; /* integers with problem dimensionality */
    int nmtt, nmty, ntsigmat, ntsigmay; /* integers with problem dimensionality */
    
    ExpData edata; /* returned rdata struct */
    UserData udata; /** user udata */
    
    udata = (UserData) user_data;
    
    /* Return rdata structure */
    edata = (ExpData) mxMalloc(sizeof *edata);
    if (edata == NULL) return(NULL);
    
    if (data_model == LW_ONEOUTPUT) {
        if ( (ny>1) | (nt>1) ) {
            mexErrMsgTxt("Data model LW_ONEOUTPUT not allowed for more than one time-point or more than one observable!");
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
            mt = mxGetPr(mxGetField(prhs[8], 0 ,"T"));
            nmty = (int) mxGetN(mxGetField(prhs[8], 0 ,"T"));
            nmtt = (int) mxGetM(mxGetField(prhs[8], 0 ,"T"));
        } else {
            mexErrMsgTxt("Field T not specified as field in data struct!");
        }
        
        if (mxGetField(prhs[8], 0 ,"Sigma_T")) {
            tsigma = mxGetPr(mxGetField(prhs[8], 0 ,"Sigma_T"));
            ntsigmay = (int) mxGetN(mxGetField(prhs[8], 0 ,"Sigma_T"));
            ntsigmat = (int) mxGetM(mxGetField(prhs[8], 0 ,"Sigma_T"));
        } else {
            mexErrMsgTxt("Field Sigma_T not specified as field in data struct!");
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
        
        if (nmtt != nmaxroot) {
            mexErrMsgTxt("Number of time-points in event matrix does not match provided nmaxroot");
        }
        
        if (ntsigmat != nmaxroot) {
            mexErrMsgTxt("Number of time-points in event-sigma matrix does not match provided nmaxroot");
        }
        
        if (nmty != nr) {
            mexErrMsgTxt("Number of events in event matrix does not match provided nr");
        }
        
        if (ntsigmay != nr) {
            mexErrMsgTxt("Number of events in event-sigma matrix does not match provided nr");
        }
        
    }


    return(edata);
}


void getRootDataFSA(int *status, int *nroots, void *ami_mem, void  *user_data, void *return_data, void *temp_data) {
    /**
     * getRootDataFSA extracts root information for forward sensitivity analysis
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
    int ir;
    
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
    if (*nroots<nmaxroot) {
        for (ir=0; ir<nr; ir++){ /* only look for roots of the rootfunction not discontinuities */
            if(rootsfound[ir] != 0) {
                rootdata[*nroots + nmaxroot*ir] = t;
                rootvaldata[*nroots + nmaxroot*ir] = 0;
                /* extract sensitivity information */
                if(sensi >= 1) {
                    if(sensi_meth == CW_FSA) {
                        *status = AMIGetSens(ami_mem, &t, sx);
                        if (*status != AMI_SUCCESS) return;
                        *status = fsroot(t,*nroots,rootSdata,x,sx,udata);
                        if (*status != AMI_SUCCESS) return;
                        if (sensi >= 2) {
                            fs2root(t,*nroots,rootS2data,x,sx,udata);
                        }
                        for (ip=0; ip<np; ip++) {
                            rootvalSdata[*nroots + nmaxroot*(ip*nr + ir)] = 0;
                            if (sensi >= 2) {
                                for (jp=0; jp<np; jp++) {
                                    rootvalS2data[*nroots + nmaxroot*((np*ip+jp)*nr + ir)] = 0;
                                }
                            }
                        }
                    }
                }
                (*nroots)++;
            }
        }
    }
    /* ROOTS FOR DISCONTINUITIES */
    for (ir=nr; ir<(nr+ndisc); ir++) {
        if(rootsfound[ir] != 0) {
            /* take care of deltas */
            /* sdeltadisc updates both x and sx, this needs to be done in combination for second order sensitivities */
            if(sensi >= 1) {
                if(sensi_meth == CW_FSA) {
                    *status = AMIGetSens(ami_mem, &t, sx);
                    if (*status != AMI_SUCCESS) return;
                    *status = sdeltadisc(t,ir-nr,x,sx,udata);
                    if (*status != AMI_SUCCESS) return;
                    *status = AMISensReInit(ami_mem, ism, sx, sdx);
                    if (*status != AMI_SUCCESS) return;
                } else {
                    *status = deltadisc(t,ir-nr,x,udata);
                    if (*status != AMI_SUCCESS) return;
                }
            } else {
                *status = deltadisc(t,ir-nr,x,udata);
                if (*status != AMI_SUCCESS) return;
            }
            /* reinitialize */
            *status = AMIReInit(ami_mem,t,x,dx);
            if (*status != AMI_SUCCESS) return;
            
            *status = AMICalcIC(ami_mem, tstart);
            if (*status != AMI_SUCCESS) return;
        }
    }
    /* if the root coincides with one of the output timepoints, we want to continue as if we reached that timepoint after accounting for the root */
}

void getRootDataASA(int *status, int *nroots, int *idisc, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data) {
    /**
     * getRootDataASA extracts root information for adjoint sensitivity analysis
     *
     * @param[out] status flag indicating success of execution @type *int
     * @param[out] nroots counter for the number of found roots @type *int
     * @param[out] idisc counter for the number of found discontinuities @type *int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] user_data pointer to the user data struct @type UserData
     * @param[out] return_data pointer to the return data struct @type ReturnData
     * @param[in] exp_data pointer to the experimental data struct @type ExpData
     * @param[out] temp_data pointer to the temporary data struct @type TempData
     * @return void
     */
    int ip;
    int ir;
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
    if (*nroots<nmaxroot) {
        for (ir=0; ir<nr; ir++){
            if(rootsfound[ir] == 1) {
                rootdata[*nroots + nmaxroot*ir] = t;
                rootvaldata[*nroots + nmaxroot*ir] = 0;
                rootidx[*nroots] = ir;
                x_tmp = NV_DATA_S(x);
                /* extract sensitivity information */
                if(!mxIsNaN(mt[ir*nmaxroot+*nroots])) {
                    if (event_model == LW_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+*nroots],2)) + 0.5*pow( ( t - mt[ir*nmaxroot+*nroots] )/tsigma[ir*nmaxroot+*nroots] , 2);
                        *chi2data += pow( ( t - mt[ir*nmaxroot+*nroots] )/tsigma[ir*nmaxroot+*nroots] , 2);
                    }
                    if (event_model == LW_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+*nroots],2)) + 0.5*pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+*nroots] , 2);
                    }
                    if (sensi>=1) {
                        x_tmp = NV_DATA_S(x);
                        *status = fdtdp(t,dtdp,x_tmp,udata);
                        if (*status != AMI_SUCCESS) return;
                        *status = fdtdx(t,dtdx,x_tmp,udata);
                        if (*status != AMI_SUCCESS) return;
                        *status = fdrvaldp(t,drvaldp,x_tmp,udata);
                        if (*status != AMI_SUCCESS) return;
                        *status = fdrvaldx(t,drvaldx,x_tmp,udata);
                        if (*status != AMI_SUCCESS) return;
                        if (mxIsNaN(tsigma[ir*nmaxroot + *nroots])) {
                            *status = fsigma_t(t,sigma_t,udata);
                            if (*status != AMI_SUCCESS) return;
                            *status = fdsigma_tdp(t,dsigma_tdp,udata);
                            if (*status != AMI_SUCCESS) return;
                        } else {
                            for (ip=0; ip<np; ip++) {
                                dsigma_ydp[ip*nr+ir] = 0;
                            }
                            sigma_t[ir] = tsigma[ir*nmaxroot + *nroots];
                        }
                        
                        for (ip=0; ip<np; ip++) {
                            if(event_model == LW_NORMAL) {
                                drdp[ip] += dsigma_tdp[ip*nr+ir]/sigma_t[ir] + ( dtdp[ip*nr+ir]* ( t - mt[ir*nmaxroot+*nroots] ) )/pow( sigma_t[ir] , 2) - dsigma_tdp[ip*nr+ir]*pow( ( t - mt[ir*nmaxroot+*nroots] ),2)/pow( sigma_t[ir] , 3);
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(event_model  == LW_NORMAL) {
                                drdx[*nroots+ix*nmaxroot] += ( dtdx[ix*nr+ir] * ( t - mt[ir*nmaxroot+*nroots] ) )/pow( sigma_t[ir] , 2);
                            }
                        }
                        for (ip=0; ip<np; ip++) {
                            if(event_model  == LW_NORMAL) {
                                drdp[ip] += dsigma_tdp[ip*nr+ir]/sigma_t[ir] + ( drvaldp[ip*nr+ir]* rootvaltmp[ir] )/pow( sigma_t[ir] , 2) - dsigma_tdp[ip*nr+ir]*pow(rootvaltmp[ir],2)/pow( sigma_t[ir] , 3);
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(event_model  == LW_NORMAL) {
                                drdx[*nroots+ix*nmaxroot] += ( drvaldx[ix*nr+ir] * rootvaltmp[ir] )/pow( sigma_t[ir] , 2);
                            }
                        }
                    }
                }
                (*nroots)++;
            }
        }
    }
    /* ROOTS FOR DISCONTINUITIES */
    for (ir=nr; ir<(nr+ndisc); ir++) {
        if(rootsfound[ir] != 0) {
            *status = deltadisc(t,ir-nr,x,udata);
            if (*status != AMI_SUCCESS) return;
            *status = AMIReInit(ami_mem,t,x,dx);
            if (*status != AMI_SUCCESS) return;
            
            *status = AMICalcIC(ami_mem, tstart);
            if (*status != AMI_SUCCESS) return;
            
            if(*idisc<nmaxdisc) {
                discs[*idisc] = t;
                irdiscs[*idisc] = ir-nr;
                (*idisc)++;
            }
        }
    }
}

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
#endif /* amici_c */