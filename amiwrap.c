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
    
    int ip, ix, it, ig; /* integers for indexing in loops */
    int ncheck; /* the number of (internal) checkpoints stored so far */
    
    void *ami_mem; /* pointer to cvodes memory block */
    UserData udata; /* user data */
    ReturnData rdata; /* return data */
    ExpData edata; /* experimental data */
    TempData tdata; /* temporary data */
    int status; /* general status flag */
    double *pstatus; /* return status flag */
    
    realtype tlastroot; /* storage for last found root */
    int iroot = 0;
    double tnext;
    
    booleantype silent;
    
    pstatus = mxMalloc(sizeof(double));
    
    udata = setupUserData(prhs);
    if (udata == NULL) goto freturn;
    
    /* options */
    if (!prhs[3]) {
        mexErrMsgIdAndTxt("AMICI:mex:options","No options provided!");
    }
    
    tdata = (TempData) mxMalloc(sizeof *tdata);
    if (tdata == NULL) goto freturn;
    
    if (nx>0) {
        ami_mem = setupAMI(&status, udata, tdata);
        if (ami_mem == NULL) goto freturn;
    }

    rdata = setupReturnData(plhs, udata, pstatus);
    if (rdata == NULL) goto freturn;
    
    if (nx>0) {
        edata = setupExpData(prhs, udata);
        if (edata == NULL) goto freturn;
    }
    
    status = workForwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot);
    if(status)
        goto freturn;

    booleantype setupBdone = false;
    status = workBackwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot, &setupBdone);

freturn:
    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);
    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, *pstatus);
    *pstatus = (double) status;
}
