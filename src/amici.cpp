/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <memory>
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

#include <include/amici.h> /* amici functions */
#include <include/amici_misc.h>
#include <include/amici_exception.h>
#include <include/symbolic_functions.h>

#include <sundials/sundials_types.h> //realtype
#include <cvodes/cvodes.h> //return codes

#include <type_traits>

// ensure definitions are in sync
static_assert(AMICI_SUCCESS == CV_SUCCESS, "AMICI_SUCCESS != CV_SUCCESS");
static_assert(AMICI_DATA_RETURN == CV_TSTOP_RETURN,
              "AMICI_DATA_RETURN != CV_TSTOP_RETURN");
static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN,
              "AMICI_ROOT_RETURN != CV_ROOT_RETURN");
static_assert(AMICI_ILL_INPUT == CV_ILL_INPUT,
              "AMICI_ILL_INPUT != CV_ILL_INPUT");
static_assert(AMICI_NORMAL == CV_NORMAL, "AMICI_NORMAL != CV_NORMAL");
static_assert(AMICI_ONE_STEP == CV_ONE_STEP, "AMICI_ONE_STEP != CV_ONE_STEP");
static_assert(std::is_same<amici::realtype, realtype>::value, "Definition of realtype does not match");

namespace amici {

/** errMsgIdAndTxt is a function pointer for printErrMsgIdAndTxt  */
msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
/** warnMsgIdAndTxt is a function pointer for printWarnMsgIdAndTxt  */
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;
    

/*!
 * runAmiciSimulation is the core integration routine. It initializes the solver
 * and runs the forward and backward problem.
 *
 * @param[in] solver Solver instance
 * @param[in] edata pointer to experimental data object @type ExpData
 * @param[in] rdata pointer to return data object @type ReturnData
 * @param[in] model model specification object @type Model
 */
void runAmiciSimulation(Solver &solver, const ExpData *edata,
                       ReturnData *rdata, Model &model) {
    if (!rdata)
        throw SetupFailure("rdata was not allocated!");

    if (model.nx <= 0) {
        return;
    }

    auto fwd = std::unique_ptr<ForwardProblem>(new ForwardProblem(rdata,edata,&model,&solver));
    fwd->workForwardProblem();

    auto bwd = std::unique_ptr<BackwardProblem>(new BackwardProblem(fwd.get()));
    bwd->workBackwardProblem();
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
