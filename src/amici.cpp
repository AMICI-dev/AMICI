/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
#define _USE_MATH_DEFINES /* MS definition of PI and other constants */
#include <cmath>
#ifndef M_PI /* define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/backwardproblem.h"
#include "include/forwardproblem.h"
#include "include/rdata.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <include/amici.h> /* amici functions */
#include <include/amici_misc.h>
#include <include/symbolic_functions.h>

msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata,
                       Model *model) {
    if (!udata)
        return AMICI_ERROR_UDATA;
    if (!rdata)
        return AMICI_ERROR_RDATA;

    int status = AMICI_SUCCESS;

    if (model->nx <= 0) {
        return AMICI_ERROR_NOTHINGTODO;
    }

    TempData tdata = TempData(udata, model, rdata);

    if (status == AMICI_SUCCESS)
        status = ForwardProblem::workForwardProblem(udata, &tdata, rdata, edata,
                                                    model);
    if (status == AMICI_SUCCESS)
        status =
            BackwardProblem::workBackwardProblem(udata, &tdata, rdata, model);

    if (status == AMICI_SUCCESS)
        status = rdata->applyChainRuleFactorToSimulationResults(udata, tdata.p);

    if (status < AMICI_SUCCESS)
        rdata->invalidate();

    return status;
}

void printErrMsgIdAndTxt(const char *identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

void printWarnMsgIdAndTxt(const char *identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
