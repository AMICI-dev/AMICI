/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
/** MS definition of PI and other constants */
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI 
/** define PI if we still have no definition */
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

/** errMsgIdAndTxt is a function pointer for printErrMsgIdAndTxt  */
msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
/** warnMsgIdAndTxt is a function pointer for printWarnMsgIdAndTxt  */
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;

/*!
 * runAmiciSimulation is the core integration routine. It initializes the solver and temporary storage in tdata and 
 * runs the forward and backward problem.
 *
 * @param[in] udata pointer to user data object @type UserData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @param[in] model pointer to model specification object @type Model
 * @return status status flag indicating (un)successful execution @type int
 */
int runAmiciSimulation(const UserData *udata, const ExpData *edata, ReturnData *rdata,
                       Model *model) {
    if (!udata || udata->nx != model->nx || udata->np != model->np || udata->nk != model->nk)
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

/*!
 * printErrMsgIdAndTxt prints a specified error message associated to the specified identifier
 *
 * @param[in] identifier error identifier @type char
 * @param[in] msg error message @type char
 * @return void
 */
void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Error] %s: %s\n", identifier, msg);
}

/*!
 * printErrMsgIdAndTxt prints a specified warning message associated to the specified identifier
 *
 * @param[in] identifier warning identifier @type char
 * @param[in] msg warning message @type char
 * @return void
 */
void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...) {
    printf("[Warning] %s: %s\n", identifier, msg);
}
