/**
 * @file   amiwrap.cpp
 * @brief  core routines for mex interface
 *
 * This file defines the fuction mexFunction which is executed upon calling the mex file from matlab
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif
#include <mex.h>
#include "wrapfunctions.h" /* user functions */
#include <include/amici_interface_matlab.h> /* amici functions */


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
    int status = 0;

    UserData udata = userDataFromMatlabCall(prhs, &status);
    ReturnDataMatlab rdata(&udata);

    plhs[0] = rdata.mxsol;

    if(status == 0 && udata.nx > 0) {
        ExpData edata = expDataFromMatlabCall(prhs, &udata, &status);
        if (status == 0)
            runAmiciSimulation(&udata, &edata, &rdata, &status);
    }

    *rdata.status = (double) status;
}
