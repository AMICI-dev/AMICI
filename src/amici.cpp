/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include "amici/amici.h"

#include "amici/backwardproblem.h"
#include "amici/forwardproblem.h"
#include "amici/misc.h"

#include <sundials/sundials_types.h> //realtype
#include <cvodes/cvodes.h> //return codes

#include <type_traits>
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
 * @param solver Solver instance
 * @param edata pointer to experimental data object
 * @param model model specification object
 * @return rdata pointer to return data object
 */
std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver, const ExpData *edata, Model &model) {
    
    auto rdata = std::unique_ptr<ReturnData>(new ReturnData(solver,&model));
    
    if (model.nx <= 0) {
        return rdata;
    }

    try{
        auto fwd = std::unique_ptr<ForwardProblem>(new ForwardProblem(rdata.get(),edata,&model,&solver));
        fwd->workForwardProblem();

        auto bwd = std::unique_ptr<BackwardProblem>(new BackwardProblem(fwd.get()));
        bwd->workBackwardProblem();
    
        rdata->status = AMICI_SUCCESS;
    } catch (amici::IntegrationFailure& ex) {
        rdata->invalidate(ex.time);
        rdata->status = ex.error_code;
        amici::warnMsgIdAndTxt("AMICI:mex:simulation","AMICI simulation failed at t = %f:\n%s\nError occured in:\n%s",ex.time,ex.what(),ex.getBacktrace());
    } catch (amici::AmiException& ex) {
        rdata->invalidate(model.t0());
        rdata->status = AMICI_ERROR;
        amici::warnMsgIdAndTxt("AMICI:mex:simulation","AMICI simulation failed:\n%s\nError occured in:\n%s",ex.what(),ex.getBacktrace());
    } catch (std::exception& ex) {
        rdata->invalidate(model.t0());
        rdata->status = AMICI_ERROR;
        amici::warnMsgIdAndTxt("AMICI:mex:simulation","AMICI simulation failed:\n%s",ex.what());
    } catch (...) {
        rdata->invalidate(model.t0());
        rdata->status = AMICI_ERROR;
        amici::warnMsgIdAndTxt("AMICI:mex", "Unknown internal error occured");
    }
    rdata->applyChainRuleFactorToSimulationResults(&model);
    return rdata;
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
