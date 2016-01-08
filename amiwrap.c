/**
 * @file   amiwrap.c
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the mex file from matlab
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
#include <src/amici.c>

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
    
    int ip, ix, jx, iy, it, ir, ie; /* integers for indexing in loops */
    int jp; /* integers for indexing in loops */
    int ncheck; /* the number of (internal) checkpoints stored so far */
    
    void *ami_mem; /* pointer to cvodes memory block */
    UserData udata; /* user data */
    ReturnData rdata; /* return data */
    ExpData edata; /* experimental data */
    TempData tdata; /* temporary data */
    int status; /* general status flag */
    double *pstatus; /* return status flag */
    int cv_status; /* status flag returned by integration method */
    
    realtype tlastroot; /* storage for last found root */
    int iroot;
    double tnext;
    
    bool rootflag, discflag;
    
    bool setupBdone;
    
    setupBdone = false;
    
    udata = setupUserData(prhs);
    if (udata == NULL) goto freturn;
    
    /* solution struct */
    
    if (!prhs[0]) {
        mexErrMsgTxt("No solution struct provided!");
    }
    
    /* options */
    if (!prhs[4]) {
        mexErrMsgTxt("No options provided!");
    }
    
    tdata = (TempData) mxMalloc(sizeof *tdata);
    if (tdata == NULL) goto freturn;
    
    ami_mem = setupAMI(&status, udata, tdata);
    if (ami_mem == NULL) goto freturn;
    
    rdata = setupReturnData(prhs, udata);
    if (rdata == NULL) goto freturn;
    
    edata = setupExpData(prhs, udata);
    if (edata == NULL) goto freturn;
    
    cv_status = 0;
    
    /*******************/
    /* FORWARD PROBLEM */
    /*******************/
    
    ncheck = 0; /* the number of (internal) checkpoints stored so far */
    
    iroot = 0;
    
    tlastroot = 0;
    /* loop over timepoints */
    for (it=0; it < nt; it++) {
        status = CVodeSetStopTime(ami_mem, ts[it]);
        if (cv_status == 0) {
            /* only integrate if no errors occured */
            if(ts[it] > tstart) {
                while (t<ts[it]) {
                    if(sensi_meth == AMI_ASA && sensi >= 1) {
                        cv_status = AMISolveF(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL, &ncheck);
                    } else {
                        cv_status = AMISolve(ami_mem, RCONST(ts[it]), x, dx, &t, AMI_NORMAL);
                    }
                    if (cv_status == -22) {
                        /* clustering of roots => turn off rootfinding */
                        AMIRootInit(ami_mem, 0, NULL);
                        cv_status = 0;
                    }
                    /* integration error occured */
                    if (cv_status<0) {
                        status = cv_status;
                        goto freturn;
                    }
                    
                    if (cv_status==AMI_ROOT_RETURN) {
                        discs[iroot] = t;
                        
                        handleEvent(&status, iroot, &tlastroot, ami_mem, udata, rdata, edata, tdata);
                        if (status != AMI_SUCCESS) goto freturn;
                        
                        if (t==ts[it]) {
                            cv_status = 0;
                        }
                        
                        iroot++;
                    }
                }
            }
            
            handleDataPoint(&status, it, ami_mem, udata, rdata, edata, tdata);
            if (cv_status==AMI_DATA_RETURN) {
                cv_status = 0;
            }
            if (status != AMI_SUCCESS) goto freturn;

            
        } else {
            for(ix=0; ix < nx; ix++) xdata[ix*nt+it] = mxGetNaN();
        }
    }

    /* fill events */
    if (ne>0) {
        fillEventOutput(&status, ami_mem, udata, rdata, edata, tdata);
    }
    
    /********************/
    /* BACKWARD PROBLEM */
    /********************/
    
    if (sensi >= 1) {
        if(sensi_meth == AMI_ASA) {
            if(cv_status == 0) {
                setupAMIB(&status, ami_mem, udata, tdata);
                setupBdone = true;
                
                it = nt-2;
                iroot--;
                while (it>=0 || iroot>=0) {
                    
                    /* check if next timepoint is a discontinuity or a data-point */
                    tnext = getTnext(discs, iroot, ts, it, udata);
                    
                    if (tnext<t) {
                        cv_status = AMISolveB(ami_mem, tnext, AMI_NORMAL);
                        
                        /* get states only if we actually integrated */
                        status = AMIGetB(ami_mem, which, &t, xB, dxB);
                        if (status != AMI_SUCCESS) goto freturn;
                        status = AMIGetQuadB(ami_mem, which, &t, xQB);
                        if (status != AMI_SUCCESS) goto freturn;
                    }
                    
                    /* handle discontinuity */
                    
                    if(ne>0){
                        if (tnext == discs[iroot]) {
                            handleEventB(&status, iroot, ami_mem, udata, tdata);
                            iroot--;
                        }
                    }
                    
                    /* handle data-point */
                    
                    if (tnext == ts[it]) {
                        handleDataPointB(&status, it, ami_mem, udata, rdata, tdata);
                        it--;
                    }
                    
                    /* reinit states */
                    status = AMIReInitB(ami_mem, which, t, xB, dxB);
                    if (status != AMI_SUCCESS) goto freturn;
                    status = AMIQuadReInitB(ami_mem, which, xQB);
                    if (status != AMI_SUCCESS) goto freturn;
                    
                    status = AMICalcICB(ami_mem, which, t, xB, dxB);
                    if (status != AMI_SUCCESS) goto freturn;
                }
                
                /* we still need to integrate from first datapoint to tstart */
                
                if (t>tstart) {
                    if(cv_status == 0) {
                        if (nx>0) {
                            /* solve for backward problems */
                            cv_status = AMISolveB(ami_mem, tstart, AMI_NORMAL);
                            
                            status = AMIGetQuadB(ami_mem, which, &t, xQB);
                            if (status != AMI_SUCCESS) goto freturn;
                            status = AMIGetB(ami_mem, which, &t, xB, dxB);
                            if (status != AMI_SUCCESS) goto freturn;
                        }
                    }
                }
                
                /* evaluate initial values */
                sx = N_VCloneVectorArray_Serial(np,x);
                if (sx == NULL) goto freturn;
                
                status = fx0(x,udata);
                if (status != AMI_SUCCESS) goto freturn;
                status = fdx0(x,dx,udata);
                if (status != AMI_SUCCESS) goto freturn;
                status = fsx0(sx, x, dx, udata);
                if (status != AMI_SUCCESS) goto freturn;
                
                if(cv_status == 0) {
                    
                    xB_tmp = NV_DATA_S(xB);
                    
                    for (ip=0; ip<np; ip++) {
                        llhS0[ip] = 0.0;
                        sx_tmp = NV_DATA_S(sx[ip]);
                        for (ix = 0; ix < nx; ix++) {
                            llhS0[ip] = llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                        }
                    }
                    
                    xQB_tmp = NV_DATA_S(xQB);
                    
                    for(ip=0; ip < np; ip++) {
                        llhSdata[ip] = - llhS0[ip] - xQB_tmp[ip];
                        if (ny>0) {
                            llhSdata[ip] -= dgdp[ip];
                        }
                        if (nz>0) {
                            llhSdata[ip] -= drdp[ip];
                        }
                    }
                    
                } else {
                    for(ip=0; ip < np; ip++) {
                        llhSdata[ip] = mxGetNaN();
                    }
                }
            } else {
                for(ip=0; ip < np; ip++) {
                    llhSdata[ip] = mxGetNaN();
                }
            }
        }
    }
    
    /* evaluate likelihood */
    
    *llhdata = - g - r;
    
    status = cv_status;
    
    goto freturn;
    
freturn:
    /* Free memory */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(dx);
        N_VDestroy_Serial(xdot);
        N_VDestroy_Serial(x_old);
        N_VDestroy_Serial(dx_old);
        N_VDestroy_Serial(xdot_old);
        DestroyMat(Jtmp);
        if (ne>0) {
            if(ami_mem) mxFree(rootsfound);
            if(ami_mem) mxFree(rootvals);
            if(ami_mem) mxFree(rootidx);
            if(ami_mem)    mxFree(sigma_z);
            if(ami_mem)    mxFree(nroots);
            if(ami_mem)    mxFree(discs);
            if(ami_mem)    mxFree(h);
            
            if(ami_mem)    mxFree(deltax);
            if(ami_mem)    mxFree(deltasx);
            if(ami_mem)    mxFree(deltaxB);
            if(ami_mem)    mxFree(deltaqB);
        }
        
        if(ny>0) {
            if(sigma_y)    mxFree(sigma_y);
        }
        if (sensi >= 1) {
            if (sensi_meth == AMI_FSA) {
                N_VDestroyVectorArray_Serial(sx,np);
            }
            if (sensi_meth == AMI_ASA) {
                if(cv_status == 0) {
                    N_VDestroyVectorArray_Serial(sx,np);
                }
            }

            if (sensi_meth == AMI_FSA) {
                N_VDestroyVectorArray_Serial(sdx, np);
            }
            if (sensi_meth == AMI_ASA) {
                if(ami_mem)    mxFree(dydx);
                if(ami_mem)    mxFree(dydp);
                if(ami_mem)    mxFree(dgdp);
                if(ami_mem)    mxFree(dgdx);
                if(ami_mem)    mxFree(drdp);
                if(ami_mem)    mxFree(drdx);
                if (ne>0) {
                    if(ami_mem)    mxFree(dzdp);
                    if(ami_mem)    mxFree(dzdx);
                }
                if(ami_mem)     mxFree(llhS0);
                if(ami_mem)    mxFree(dsigma_ydp);
                if (ne>0) {
                    if(ami_mem)    mxFree(dsigma_zdp);
                }
                if(setupBdone)      N_VDestroy_Serial(dxB);
                if(setupBdone)      N_VDestroy_Serial(xB);
                if(setupBdone)     N_VDestroy_Serial(xB_old);
                if(setupBdone)      N_VDestroy_Serial(xQB);
                if(setupBdone)      N_VDestroy_Serial(xQB_old);
            }
            
        }
        if(ami_mem)     N_VDestroy_Serial(id);
        if(ami_mem)     AMIFree(&ami_mem);
    }
    
    if(udata)   mxFree(plist);
    if (sensi >= 1) {
        if(udata)   mxFree(tmp_dxdotdp);
    }
    
    if(udata)    mxFree(udata);
    if(tdata)    mxFree(tdata);
    
    if(mxGetField(prhs[0], 0 ,"status")) { pstatus = mxGetPr(mxGetField(prhs[0], 0 ,"status")); } else { mexErrMsgTxt("Parameter status not specified as field in solution struct!"); }
    *pstatus = (double) status;
    
    return;
}