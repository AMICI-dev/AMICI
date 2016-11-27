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
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>
#include <include/udata_accessors.h>
#include <include/tdata_accessors.h>
#include <include/edata_accessors.h>
#include <include/rdata_accessors.h>

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
    
    void *ami_mem = nullptr; /* pointer to cvodes memory block */
    UserData *udata = nullptr; /* user data */
    ReturnData *rdata = nullptr; /* return data */
    ExpData *edata = nullptr; /* experimental data */
    TempData *tdata = nullptr; /* temporary data */
    int status = 0; /* general status flag */
    double *pstatus; /* return status flag */

    int iroot = 0;
    booleantype setupBdone = false;

    pstatus = (double *) mxMalloc(sizeof(double));
    
    udata = setupUserData(prhs);
    if (udata == NULL) {
        /* goto freturn will fail here as freeXXXXData routines will fail*/
        *pstatus = -98;
        return;
    }
    
    /* options */
    if (!prhs[3]) {
        mexErrMsgIdAndTxt("AMICI:mex:options","No options provided!");
    }
    
    tdata = new TempData();
    if (tdata == NULL) goto freturn;
    
    if (nx>0) {
        ami_mem = setupAMI(&status, udata, tdata);
        if (ami_mem == NULL){
            status = -96;
            goto freturn;
        }
    }

    rdata = setupReturnData(plhs, udata, pstatus);
    if (rdata == NULL) {
        status = -96;
        goto freturn;
    }
    
    
    if (nx>0) {
        edata = setupExpData(prhs, udata);
        if (edata == NULL) {
            status = -97;
            goto freturn;
        }
    }
    
    if(nx>0) {
        status = workForwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot);
        if(status<0) goto freturn;

        status = workBackwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot, &setupBdone);
    }

freturn:
    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);
    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, *pstatus);
    freeUserData(udata);
    delete edata;
    *pstatus = (double) status;
}
