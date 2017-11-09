/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
/** MS definition of PI and other constants */
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
/** define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/amici_exception.h"
#include "include/backwardproblem.h"
#include "include/forwardproblem.h"
#include "include/rdata.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <include/amici.h> /* amici functions */
#include <include/amici_misc.h>
#include <include/amici_exception.h>
#include <include/symbolic_functions.h>

namespace amici {

/** errMsgIdAndTxt is a function pointer for printErrMsgIdAndTxt  */
msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
/** warnMsgIdAndTxt is a function pointer for printWarnMsgIdAndTxt  */
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;
    

/*!
 * runAmiciSimulation is the core integration routine. It initializes the solver
 * and temporary storage in tdata and
 * runs the forward and backward problem.
 *
 * @param[in] udata pointer to user data object @type UserData
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @param[in] model pointer to model specification object @type Model
 */
void runAmiciSimulation(const UserData *udata, const ExpData *edata,
                       ReturnData *rdata, Model *model) {
    if (!udata || udata->nx != model->nx || udata->np != model->np ||
        udata->nk != model->nk)
        throw SetupFailure("udata was not allocated or does not agree with model!");
    if (!rdata)
        throw SetupFailure("rdata was not allocated!");

    if (model->nx <= 0) {
        return;
    }

    TempData tdata(udata, model, rdata);

    ForwardProblem::workForwardProblem(udata, &tdata, rdata, edata,
                                                    model);
    BackwardProblem::workBackwardProblem(udata, &tdata, rdata, model);

    return;
}

/*!
 * printErrMsgIdAndTxt prints a specified error message associated to the
 * specified identifier
 *
 * @param[in] identifier error identifier @type char
 * @param[in] format string with error message printf-style format
 * @param ... arguments to be formatted
 * @return void
 */
void printErrMsgIdAndTxt(const char *identifier, const char *format, ...) {
    if(identifier != NULL && *identifier != '\0')
        fprintf(stderr, "[Error] %s: ", identifier);
    else
        fprintf(stderr, "[Error] ");
    va_list argptr;
    va_start(argptr,format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);
    fprintf(stderr, "\n");
}

/*!
 * printErrMsgIdAndTxt prints a specified warning message associated to the
 * specified identifier
 *
 * @param[in] identifier warning identifier @type char
 * @param[in] format string with error message printf-style format
 * @param ... arguments to be formatted
 * @return void
 */
void printWarnMsgIdAndTxt(const char *identifier, const char *format, ...) {
    if(identifier != NULL && *identifier != '\0')
        printf("[Warning] %s: ", identifier);
    else
        printf("[Warning] ");
    va_list argptr;
    va_start(argptr,format);
    vprintf(format, argptr);
    va_end(argptr);
    printf("\n");
}

} // namespace amici
