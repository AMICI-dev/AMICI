/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
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

/** return value for successful execution */
#define AMI_SUCCESS               0

static int fsy(realtype t_, int it, realtype *sy, realtype *dydx_, realtype *dydp_, N_Vector *sx, void *user_data);
static int fsJy(realtype t_, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y_, realtype *sy, realtype *dydp_, realtype *my_, void *user_data);

void runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata, int *pstatus) {

    if(nx <= 0) {
        *pstatus = -99;
        return;
    }

    *pstatus = 0;

    TempData *tdata = new TempData();
    if (tdata == NULL) {
        *pstatus = -100;
        return;
    }

    unscaleParameters(udata);

    /* pointer to cvodes memory block */
    void *ami_mem = setupAMI(pstatus, udata, tdata);
    if (ami_mem == NULL){
        *pstatus = -96;
        delete tdata;
        return;
    }

    int iroot = 0;
    booleantype setupBdone = false;

    int problem = workForwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot);
    if(problem)
        goto freturn;

    problem = workBackwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot, &setupBdone);
    if(problem)
        goto freturn;

    applyChainRuleFactorToSimulationResults(udata, rdata, edata);

freturn:
    if(*pstatus<0){
        invalidateReturnData(udata, rdata);
    }
    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, *pstatus);
}

void invalidateReturnData(UserData* udata, ReturnData* rdata) {
    /**
     * @brief performs all necessary actions to reset return data upon integration failure
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     */
    if(llhdata)
        *llhdata = amiGetNaN();

    if(sllhdata)
        fillArray(sllhdata, nplist, amiGetNaN());

    if(s2llhdata)
        fillArray(s2llhdata, nplist*(ng-1), amiGetNaN());
}

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
        if(ne>0) deltasx = new realtype[nx*nplist]();
        if(ne>0) deltaxB = new realtype[nx]();
        if(ne>0) deltaqB = new realtype[ng*nplist]();

        if(ny>0) sigma_y = new realtype[ny]();
        if(ne>0) sigma_z = new realtype[nz]();


        /* initialise states */
        if (x == NULL) return(NULL);
        if(x0data == NULL) {
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
        if (*status != AMI_SUCCESS) return(NULL);

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

    if (sensi >= AMI_SENSI_ORDER_FIRST) {

        dydx = new realtype[ny*nx]();
        dydp = new realtype[ny*nplist]();
        dgdp = new realtype[ng*nplist*nytrue]();
        dgdx = new realtype[ng*nxtrue*nt]();
        dgdy = new realtype[nytrue*ng*ny]();
        if (ne > 0) {
            dzdp = new realtype[nz*nplist]();
            dzdx = new realtype[nz*nx]();
        }
        drdp = new realtype[ng*nplist*nztrue*nmaxevent]();
        drdx = new realtype[ng*nx*nztrue*nmaxevent]();

        dsigma_ydp = new realtype[ny*nplist]();
        if(ne>0) dsigma_zdp = new realtype[nz*nplist]();

        if (sensi_meth == AMI_SENSI_FSA) {

            if(nx>0) {

                /* allocate some more temporary storage */
                NVsx = N_VCloneVectorArray_Serial(nplist, x);
                sdx = N_VCloneVectorArray_Serial(nplist, x);
                if (NVsx == NULL) return(NULL);
                if (sdx == NULL) return(NULL);

                /* initialise sensitivities, this can either be user provided or come from the model definition */

                if(!sx0data) {
                    *status = fsx0(NVsx, x, dx, udata);
                    if (*status != AMI_SUCCESS) return(NULL);
                } else {
                    int ip;
                    for (ip=0; ip<nplist; ip++) {
                        sx_tmp = NV_DATA_S(NVsx[ip]);
                        int ix;
                        for (ix=0; ix<nx; ix++) {
                            sx_tmp[ix] = sx0data[ix + nx*ip];
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

        if (sensi_meth == AMI_SENSI_ASA) {

            if(nx>0) {
                /* Allocate space for the adjoint computation */

                which = 0;

                if(ne>0) x_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);
                if(ne>0) xdot_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);
                if(ne>0) xdot_old_disc = N_VCloneVectorArray_Serial(ne*nmaxevent, x);

                *status = AMIAdjInit(ami_mem, maxsteps, interpType);
                if (*status != AMI_SUCCESS) return(NULL);

                llhS0 = new realtype[ng*nplist]();
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
    xQB = N_VNew_Serial(ng*nplist);
    xQB_old = N_VNew_Serial(ng*nplist);

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
    memset(xQB_tmp,0,sizeof(realtype)*ng*nplist);

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
    *status = AMISetMaxNumStepsB(ami_mem, which, 10*maxsteps);
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

void getDataSensisFSA(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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

    for(ip=0; ip < nplist; ip++) {
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
        if(edata){
            if (amiIsNaN(ysigma[iy*nt+it])) {
                *status = fdsigma_ydp(t,dsigma_ydp,udata);
                if (*status != AMI_SUCCESS) return;
            } else {
                for (ip=0; ip<nplist; ip++) {
                    dsigma_ydp[ip*ny+iy] = 0;
                }
            }
            for (ip=0; ip<nplist; ip++) {
                ssigmaydata[it + nt*(ip*ny+iy)] = dsigma_ydp[ip*ny+iy];
            }
        } else {
            for (ip=0; ip<nplist; ip++) {
                ssigmaydata[it + nt*(ip*ny+iy)] = 0;
            }
        }
    }
    fsy(ts[it],it,sydata,dydx,dydp,NVsx,udata);
    if(edata) {
        fsJy(ts[it],it,sllhdata,s2llhdata,dgdy,dgdp,ydata,sigma_y,sydata,dydp,my,udata);
    }
}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void prepDataSensis(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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
    if(edata) {
        for (iy=0; iy<nytrue; iy++) {
            if (amiIsNaN(ysigma[iy*nt+it])) {
                *status = fdsigma_ydp(t,dsigma_ydp,udata);
                if (*status != AMI_SUCCESS) return;
            } else {
                for (ip=0; ip<nplist; ip++) {
                    dsigma_ydp[ip*ny+iy] = 0;
                }
            }
            for (ip=0; ip<nplist; ip++) {
                ssigmaydata[it + nt*(ip*ny+iy)] = dsigma_ydp[ip*ny+iy];
            }
        }
        fdJydp(ts[it],it,dgdp,ydata,x,dydp,my,sigma_y,dsigma_ydp,udata);


        if (sensi_meth == AMI_SENSI_ASA) {
            for(ig=0; ig<ng; ig++) {
                for(ip=0; ip < nplist; ip++) {
                    for(iy=0; iy < nytrue; iy++) {
                        if(ig==0) {
                            if (ny>0) {
                                sllhdata[ip] -= dgdp[iy + ip*nytrue];
                            }
                        } else {
                            if (ny>0) {
                                s2llhdata[(ig-1)*nplist + ip] -= dgdp[(ig*nplist + ip)*nytrue + iy];
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

void getDataOutput(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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

    if(edata) {
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
    } else {
        *status =fsigma_y(t,sigma_y,udata);
        if (*status != AMI_SUCCESS) return;
        for (iy=0; iy<nytrue; iy++) {
            sigmaydata[iy*nt+it] = sigma_y[iy];
        }
    }
    if (sensi >= AMI_SENSI_ORDER_FIRST) {
        prepDataSensis(status, it, ami_mem, udata, rdata, edata, tdata);
        if (sensi_meth == AMI_SENSI_FSA) {
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

    if(sensi >= AMI_SENSI_ORDER_SECOND) {
        *status = fs2root(t,ie,nroots,s2rzdata,x,NVsx,udata);
        if (*status != AMI_SUCCESS) return;
    }

}

/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void getEventSensisASA(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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
                    for (ip=0; ip<nplist; ip++) {
                        dsigma_zdp[iz+nz*ip] = 0;
                    }
                    sigma_z[iz] = zsigma[nroots[ie] + nmaxevent*iz];
                }
                sigmazdata[nroots[ie] + nmaxevent*iz] = sigma_z[iz];
                for (ip=0; ip<nplist; ip++) {
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

void getEventSigma(int *status, int ie, int iz, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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

void getEventObjective(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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
    if(edata) {
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

void getEventOutput(int *status, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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
                if (sensi >= AMI_SENSI_ORDER_FIRST) {
                    if(sensi_meth == AMI_SENSI_ASA) {
                        getEventSensisASA(status, ie, ami_mem, udata, rdata, edata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    } else {
                        getEventSensisFSA(status, ie, ami_mem, udata, rdata, tdata);
                        if (*status != AMI_SUCCESS) return;
                    }
                }

                if(edata) {
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

void fillEventOutput(int *status, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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

                if (sensi >= AMI_SENSI_ORDER_FIRST) {
                    if(sensi_meth == AMI_SENSI_ASA) {
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

void handleDataPoint(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata) {
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
            if( sensi_meth == AMI_SENSI_SS) {

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

void handleEvent(int *status, int *iroot, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, int seflag) {
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
        if(sensi >= AMI_SENSI_ORDER_FIRST){
            if (sensi_meth == AMI_SENSI_FSA) {
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
    if(sensi >= AMI_SENSI_ORDER_FIRST){
        /* store x and xdot to compute jump in sensitivities */
        N_VScale(1.0,x,x_old);
        if (sensi_meth == AMI_SENSI_FSA) {
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

        if (sensi_meth == AMI_SENSI_ASA) {
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

    if(sensi >= AMI_SENSI_ORDER_FIRST){
        if (sensi_meth == AMI_SENSI_FSA) {

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

    if(sensi >= AMI_SENSI_ORDER_FIRST){
        if (sensi_meth == AMI_SENSI_FSA) {
            if (seflag == 0) {
                *status = AMISensReInit(ami_mem, ism, NVsx, sdx);
                if (*status != AMI_SUCCESS) return;
            }
        }
    }
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
                for (ip=0; ip<nplist; ip++) {
                    xQB_tmp[ig*nplist+ip] += deltaqB[ig*nplist+ip];
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

            for (ip=0; ip<nplist; ip++) {
                sx_tmp = NV_DATA_S(NVsx[ip]);
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
        if (rootvals[ie]<0) {
            h[ie] = 0.0;
        } else if (rootvals[ie]==0) {
            errMsgIdAndTxt("AMICI:mex:initHeaviside","Simulation started in an event. This could lead to unexpected results, aborting simulation! Please specify an earlier simulation start via @amimodel.t0");
            *status = -10;
            return;
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


int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, int *status, void *ami_mem, int *iroot) {
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
        if(sensi_meth == AMI_SENSI_FSA && sensi >= AMI_SENSI_ORDER_FIRST) {
            *status = AMISetStopTime(ami_mem, ts[it]);
        }
        if (*status == 0) {
            /* only integrate if no errors occured */
            if(ts[it] > tstart) {
                while (t<ts[it]) {
                    if(sensi_meth == AMI_SENSI_ASA && sensi >= AMI_SENSI_ORDER_FIRST) {
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

int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, int *status, void *ami_mem, int *iroot, booleantype *setupBdone) {
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
        if (sensi >= AMI_SENSI_ORDER_FIRST) {
            if(sensi_meth == AMI_SENSI_ASA) {
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
                    NVsx = N_VCloneVectorArray_Serial(nplist,x);
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
                                for (ip=0; ip<nplist; ip++) {
                                    llhS0[ig*nplist + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(NVsx[ip]);
                                    for (ix = 0; ix < nxtrue; ix++) {
                                        llhS0[ip] = llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                                    }
                                }
                            } else {
                                for (ip=0; ip<nplist; ip++) {
                                    llhS0[ig*nplist + ip] = 0.0;
                                    sx_tmp = NV_DATA_S(NVsx[ip]);
                                    for (ix = 0; ix < nxtrue; ix++) {
                                        llhS0[ig*nplist + ip] = llhS0[ig*nplist + ip] + xB_tmp[ig*nxtrue + ix] * sx_tmp[ix] + xB_tmp[ix] * sx_tmp[ig*nxtrue + ix];
                                    }
                                }
                            }
                        }

                        xQB_tmp = NV_DATA_S(xQB);

                        for(ig=0; ig<ng; ig++) {
                            for(ip=0; ip < nplist; ip++) {
                                if (ig==0) {
                                    sllhdata[ip] -=  llhS0[ip] + xQB_tmp[ip];
                                    if (nz>0) {
                                        sllhdata[ip] -= drdp[ip];
                                    }
                                } else {
                                    s2llhdata[(ig-1)*nplist + ip] -= llhS0[ig*nplist + ip] + xQB_tmp[ig*nplist + ip];
                                    if (nz>0) {
                                        s2llhdata[(ig-1)*nplist + ip] -= drdp[ig*nplist + ip];
                                    }
                                }
                            }
                        }

                    } else {
                        int ig;
                        for(ig=0; ig<ng; ig++) {
                            for(ip=0; ip < nplist; ip++) {
                                if (ig==0) {
                                    sllhdata[ip] = amiGetNaN();
                                } else {
                                    s2llhdata[(ig-1)*nplist + ip] = amiGetNaN();
                                }
                            }
                        }
                    }
                } else {
                    int ig;
                    for(ig=0; ig<ng; ig++) {
                        for(ip=0; ip < nplist; ip++) {
                            if (ig==0) {
                                sllhdata[ip] = amiGetNaN();
                            } else {
                                s2llhdata[(ig-1)*nplist + ip] = amiGetNaN();
                            }
                        }
                    }
                }
            }
        }
    }

    /* evaluate likelihood */
    if(edata) {
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
                if(xdotdata)
                    if(xdot_tmp)
                        memcpy(xdotdata,xdot_tmp,nx*sizeof(realtype));
            }
        }
    }
    if(udata) {
        if(nx>0) {
            fJ(nx,t,0,x,dx,xdot,Jtmp,udata,NULL,NULL,NULL);
            if(Jdata)
                if(Jtmp->data)
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
        if (sensi >= AMI_SENSI_ORDER_FIRST) {
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
            if (sensi_meth == AMI_SENSI_FSA) {
                N_VDestroyVectorArray_Serial(NVsx,nplist);
            }
            if (sensi_meth == AMI_SENSI_ASA) {
                if(NVsx) {
                    N_VDestroyVectorArray_Serial(NVsx,nplist);
                }
            }

            if (sensi_meth == AMI_SENSI_FSA) {
                N_VDestroyVectorArray_Serial(sdx, nplist);
            }
            if (sensi_meth == AMI_SENSI_ASA) {

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


void unscaleParameters(UserData *udata) {
    switch(udata->am_pscale) {
    case AMI_SCALING_LOG10:
        for(int ip = 0; ip < np; ++ip) {
            p[ip] = pow(10, p[ip]);
        }
        break;
    case AMI_SCALING_LN:
        for(int ip = 0; ip < np; ++ip)
            p[ip] = exp(p[ip]);
        break;
    case AMI_SCALING_NONE:
        //this should never be reached
        break;
    }
}

void applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata, const ExpData *edata)
{
    if(udata->am_pscale == AMI_SCALING_NONE)
        return;

    // chain-rule factor: multiplier for am_p
    realtype coefficient;
    realtype *pcoefficient, *augcoefficient;

    pcoefficient = new realtype[nplist]();
    augcoefficient = new realtype[np]();

    switch(udata->am_pscale) {
    case AMI_SCALING_LOG10:
            coefficient = log(10.0);
            for(int ip = 0; ip < nplist; ++ip)
                pcoefficient[ip] = p[plist[ip]]*log(10);
            if(sensi == 2)
                if(udata->am_o2mode == AMI_O2MODE_FULL)
                    for(int ip = 0; ip < np; ++ip)
                    augcoefficient[ip] = p[ip]*log(10);
        break;
    case AMI_SCALING_LN:
            coefficient = 1.0;
            for(int ip = 0; ip < nplist; ++ip)
                pcoefficient[ip] = p[plist[ip]];
            if(sensi == 2)
                if(udata->am_o2mode == AMI_O2MODE_FULL)
                    for(int ip = 0; ip < np; ++ip)
                        augcoefficient[ip] = p[ip];
        break;
    case AMI_SCALING_NONE:
            //this should never be reached
        break;
    }

    if(sensi >= AMI_SENSI_ORDER_FIRST) {
        // recover first order sensitivies from states for adjoint sensitivity analysis
        if(sensi == AMI_SENSI_ORDER_SECOND){
            if(sensi_meth == AMI_SENSI_ASA){
                if(rdata->am_xdata)
                    if(rdata->am_sxdata)
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int ix = 0; ix < nxtrue; ++ix)
                                for(int it = 0; it < nt; ++it)
                                    sxdata[(ip*nxtrue + ix)*nt + it] = xdata[(nxtrue + ip*nxtrue + ix)*nt + it];

                if(rdata->am_ydata)
                    if(rdata->am_sydata)
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int iy = 0; iy < nytrue; ++iy)
                                for(int it = 0; it < nt; ++it)
                                    sydata[(ip*nytrue + iy)*nt + it] = ydata[(nytrue + ip*nytrue + iy)*nt + it];

                if(rdata->am_zdata)
                    if(rdata->am_szdata)
                        for(int ip = 0; ip < nplist; ++ip)
                            for(int iz = 0; iz < nztrue; ++iz)
                                for(int it = 0; it < nt; ++it)
                                    sydata[(ip*nztrue + iz)*nt + it] = zdata[(nztrue + ip*nztrue + iz)*nt + it];

            }
        }

        if(edata) {
            if(rdata->am_sllhdata)
                for(int ip = 0; ip < nplist; ++ip)
                    sllhdata[ip] *= pcoefficient[ip];
        }

#define chainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if(rdata->am_s ## QUANT ## data) \
    for(int ip = 0; ip < nplist; ++ip) \
        for(int IND1 = 0; IND1 < N1T; ++IND1) \
            for(int IND2 = 0; IND2 < N2; ++IND2){ \
                s ## QUANT ## data[(ip * N1 + IND1) * N2 + IND2] *= pcoefficient[ip];} \

        chainRule(x,ix,nxtrue,nx,it,nt)
        chainRule(y,iy,nytrue,ny,it,nt)
        chainRule(sigmay,iy,nytrue,ny,it,nt)
        chainRule(z,iz,nztrue,nz,ie,nmaxevent)
        chainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        chainRule(rz,iz,nztrue,nz,ie,nmaxevent)
    }
    if(sensi_meth == AMI_SENSI_SS) {
        if(rdata->am_dxdotdpdata)
            for(int ip = 0; ip < nplist; ++ip)
                for(int ix = 0; ix < nx; ++ix)
                    dxdotdpdata[ip*nxtrue + ix] *= pcoefficient[ip];

        if(rdata->am_dydpdata)
            for(int ip = 0; ip < nplist; ++ip)
                for(int iy = 0; iy < ny; ++iy)
                    dydpdata[ip*nxtrue + iy] *= pcoefficient[ip];
    }
    if(udata->am_o2mode == AMI_O2MODE_FULL) { //full
        if(edata){
            if(rdata->am_s2llhdata) {
                if(rdata->am_sllhdata) {
                    for(int ip = 0; ip < nplist; ++ip) {
                        for(int ig = 1; ig < ng; ++ig) {
                            s2llhdata[ip*nplist+(ig-1)] *= pcoefficient[ip]*augcoefficient[ig-1];
                            if(plist[ip] == ig-1)
                                s2llhdata[ip*nplist+(ig-1)] += sllhdata[ip]*coefficient;
                        }
                    }
                }
            }
        }

#define s2ChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if(rdata->am_s ## QUANT ## data) \
    for(int ip = 0; ip < nplist; ++ip) \
        for(int ig = 1; ig < ng; ++ig) \
            for(int IND1 = 0; IND1 < N1T; ++IND1) \
                for(int IND2 = 0; IND2 < N2; ++IND2){ \
                    s ## QUANT ## data[(ip*N1 + ig*N1T + IND1)*N2 + IND2] *= pcoefficient[ip]*augcoefficient[ig-1]; \
                    if(plist[ip]==ig-1) \
                        s  ## QUANT ## data[(ip*N1 + ig*N1T + IND1)*N2 + IND2] += s ## QUANT ## data[(ip*N1 + IND1)*N2 + IND2]*coefficient;}

        s2ChainRule(x,ix,nxtrue,nx,it,nt)
        s2ChainRule(y,iy,nytrue,ny,it,nt)
        s2ChainRule(sigmay,iy,nytrue,ny,it,nt)
        s2ChainRule(z,iz,nztrue,nz,ie,nmaxevent)
        s2ChainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        s2ChainRule(rz,iz,nztrue,nz,ie,nmaxevent)
    }

    if(udata->am_o2mode == AMI_O2MODE_DIR) { //directional
        if(rdata->am_s2llhdata) {
            if(rdata->am_sllhdata) {
                for(int ip = 0; ip < nplist; ++ip) {
                    s2llhdata[ip] *= pcoefficient[ip];
                    s2llhdata[ip] += udata->am_k[nk-nplist+ip]*sllhdata[ip]/p[plist[ip]];
                }
            }
        }

#define s2vecChainRule(QUANT,IND1,N1T,N1,IND2,N2) \
if(rdata->am_s ## QUANT ## data) \
    for(int ip = 0; ip < nplist; ++ip) \
            for(int IND1 = 0; IND1 < N1T; ++IND1) \
                for(int IND2 = 0; IND2 < N2; ++IND2){ \
                    s ## QUANT ## data[(ip*N1 + N1T + IND1)*N2 + IND2] *= pcoefficient[ip]; \
                    s ## QUANT ## data[(ip*N1 + N1T + IND1)*N2 + IND2] += udata->am_k[nk-nplist+ip]*s ## QUANT ## data[(ip*N1 + IND1)*N2 + IND2]/p[plist[ip]];}

        s2vecChainRule(x,ix,nxtrue,nx,it,nt)
        s2vecChainRule(y,iy,nytrue,ny,it,nt)
        s2vecChainRule(sigmay,iy,nytrue,ny,it,nt)
        s2vecChainRule(z,iz,nztrue,nz,ie,nmaxevent)
        s2vecChainRule(sigmaz,iz,nztrue,nz,ie,nmaxevent)
        s2vecChainRule(rz,iz,nztrue,nz,ie,nmaxevent)
    }

    delete[] pcoefficient;
    delete[] augcoefficient;
}

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
    if (sensi >= AMI_SENSI_ORDER_FIRST) {
        /* initialise temporary dxdotdp storage */
        tmp_dxdotdp = new realtype[nx*nplist]();
    }
    if (ne>0) {
        /* initialise temporary stau storage */
        stau_tmp = new realtype[nplist]();
    }


    w_tmp = new realtype[nw]();
    dwdx_tmp = new realtype[ndwdx]();
    dwdp_tmp = new realtype[ndwdp]();
}

int fsy(realtype t_, int it, realtype *sy, realtype *dydx_, realtype *dydp_, N_Vector *sx, void *user_data){
    // Compute sy = dydx * sx + dydp

    int status = 0;
    UserData *udata = (UserData*) user_data;

    for(int ip = 0; ip < nplist; ++ip) {
        for(int iy = 0; iy < ny; ++iy)
            // copy dydp to sy
            sy[ip * nt * ny + iy * nt + it] = dydp_[iy + ip * ny];

        realtype *sxTmp = N_VGetArrayPointer(sx[ip]);

        // compute sy = 1.0*dydx*sx + 1.0*sy
        cblas_dgemv(CblasColMajor, CblasNoTrans, ny, nx,
                    1.0, dydx_, ny, sxTmp, 1,
                    1.0, &sy[ip * nt * ny + it], nt);
    }

    return status;
}

int fsJy(realtype t_, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y_, realtype *sy, realtype *dydp_, realtype *my_, void *user_data) {
    int status = 0;
    UserData *udata = (UserData*) user_data;

    // Compute sy-dydp for current 'it'
    // dydp         ny x nplist
    // sy           nt x ny x nṕlist
    // dydp part needs to be substracted as it is already contained in dJydp
    // we only need to account for sensitivities here
    realtype *diff = new realtype[ny * nplist];
    for(int iy = 0; iy < ny; ++iy)
        for(int ip = 0; ip < nplist; ++ip)
            diff[iy + ip * ny] = sy[ip * nt * ny + iy * nt + it] - dydp_[iy + ip * ny];

    // sJy          nplist x ng
    // dJydp=dgdp   nytrue x nplist x ng
    // dJydy=dgdy   nytrue x ng x ny

    realtype *dJydyTmp = new realtype[ng * ny];
    realtype *multResult = new realtype[nplist * ng];

    for(int iyt = 0; iyt < nytrue; ++iyt) {
        if(amiIsNaN(my_[nt * iyt + it]))
            continue;

        // copy current (iyt) dJydy slice
        // dJydyTmp     ng x ny
        for(int ig = 0; ig < ng; ++ig)
            for(int iy = 0; iy < ny; ++iy)
                dJydyTmp[ig + iy * ng] = dJydy[iyt + ig * nytrue + iy * nytrue * ng];

        // compute multResult = (dJydyTmp * diff)' + dJydp == diff' * dJydyTmp' + dJydp
        // copy dJydp slice (iyt) to result
        for(int ip = 0; ip < nplist; ++ip)
            for(int ig = 0; ig < ng; ++ig)
                multResult[ip + np * ig] = dJydp[iyt + ip * nytrue + ig * nytrue * nplist];

        // C := alpha*op(A)*op(B) + beta*C,
        cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans,
                    nplist, ng, ny,
                    1.0, diff, ny,
                    dJydyTmp, ng,
                    1.0, multResult, nplist);


        // sJy += multResult
        for(int ig = 0; ig < ng; ++ig) {
            if(ig == 0)
                for(int ip = 0; ip < nplist; ++ip)
                    sJy[ip] -= multResult[ip];
            else
                for(int ip = 0; ip < nplist; ++ip)
                    s2Jy[ip + nplist * (ig - 1)] -= multResult[ip+ nplist * ig];
        }


    }
    delete[] dJydyTmp;
    delete[] multResult;
    delete[] diff;

    return(status);
}

