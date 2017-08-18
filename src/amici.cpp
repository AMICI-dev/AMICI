/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <include/amici.h> /* amici functions */
#include <include/symbolic_functions.h>
#include <include/amici_misc.h>
#include "include/amici_solver.h"
#include "include/amici_model.h"
#include "include/forwardproblem.h"
#include "include/backwardproblem.h"

msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata, Model *model, Solver *solver) {
    if(!udata) return AMICI_ERROR_UDATA;
    if(!rdata) return AMICI_ERROR_RDATA;
    
    int status = AMICI_SUCCESS;
    
    if (model->nx <= 0) {
        return AMICI_ERROR_NOTHINGTODO;
    }
    
    TempData *tdata = new TempData(udata, model);
    
    // unscale parameters but keep original
    double *originalParams = NULL;
    if(model->pscale != AMICI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * model->np);
        memcpy(originalParams, udata->p, sizeof(double) * model->np);
    }
    status = udata->unscaleParameters(model);

    if (status == AMICI_SUCCESS)
        status = solver->setupAMI(udata, tdata, model);

    if (status != AMICI_SUCCESS)
        goto freturn;

    if (status == AMICI_SUCCESS) status = ForwardProblem::workForwardProblem(udata, tdata, rdata, edata, solver, model);
    if (status == AMICI_SUCCESS) status = BackwardProblem::workBackwardProblem(udata, tdata, rdata, solver, model);
    
    if (status == AMICI_SUCCESS) status = rdata->applyChainRuleFactorToSimulationResults(udata);
    if (status < AMICI_SUCCESS) rdata->invalidate();
    
    
freturn:
    // reset to original parameters
    if(originalParams) {
        memcpy(udata->p, originalParams, sizeof(double) * model->np);
        free(originalParams);
    }

    delete tdata;

    return status;
}


/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------- */

void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
