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
static_assert(std::is_same<amici::realtype, realtype>::value,
              "Definition of realtype does not match");


namespace amici {

msgIdAndTxtFp errMsgIdAndTxt = &printErrMsgIdAndTxt;
msgIdAndTxtFp warnMsgIdAndTxt = &printWarnMsgIdAndTxt;


std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver, const ExpData *edata, Model &model, bool rethrow) {
    std::unique_ptr<ReturnData> rdata;

    /* Applies condition-specific model settings and restores them when going
     * out of scope */
    ConditionContext conditionContext(&model, edata);

    try{
        rdata = std::unique_ptr<ReturnData>(new ReturnData(solver, model));

        if (model.nx_solver <= 0) {
            return rdata;
        }

        auto fwd = std::unique_ptr<ForwardProblem>(new ForwardProblem(rdata.get(),edata,&model,&solver));
        fwd->workForwardProblem();

        auto bwd = std::unique_ptr<BackwardProblem>(new BackwardProblem(fwd.get()));
        bwd->workBackwardProblem();

        rdata->status = AMICI_SUCCESS;
    } catch (amici::IntegrationFailure const& ex) {
        rdata->invalidate(ex.time);
        rdata->status = ex.error_code;
        if(rethrow) throw;
        amici::warnMsgIdAndTxt("AMICI:mex:simulation","AMICI forward simulation failed at t = %f:\n%s\n",ex.time,ex.what());
    } catch (amici::IntegrationFailureB const& ex) {
        rdata->invalidateSLLH();
        rdata->status = ex.error_code;
        if(rethrow) throw;
        amici::warnMsgIdAndTxt(
                    "AMICI:mex:simulation",
                    "AMICI backward simulation failed when trying to solve until t = %f"
                    " (see message above):\n%s\n",
                    ex.time, ex.what());
    } catch (amici::AmiException const& ex) {
        rdata->invalidate(model.t0());
        rdata->status = AMICI_ERROR;
        if(rethrow) throw;
        amici::warnMsgIdAndTxt("AMICI:mex:simulation","AMICI simulation failed:\n%s\nError occured in:\n%s",ex.what(),ex.getBacktrace());
    }

    rdata->applyChainRuleFactorToSimulationResults(&model);

    return rdata;
}

void printErrMsgIdAndTxt(const char *identifier, const char *format, ...) {
    if(identifier && *identifier != '\0')
        fprintf(stderr, "[Error] %s: ", identifier);
    else
        fprintf(stderr, "[Error] ");
    va_list argptr;
    va_start(argptr,format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);
    fprintf(stderr, "\n");
}

void printWarnMsgIdAndTxt(const char *identifier, const char *format, ...) {
    if(identifier && *identifier != '\0')
        printf("[Warning] %s: ", identifier);
    else
        printf("[Warning] ");
    va_list argptr;
    va_start(argptr,format);
    vprintf(format, argptr);
    va_end(argptr);
    printf("\n");
}

std::vector<std::unique_ptr<ReturnData> > runAmiciSimulations(const Solver &solver,
                                                              const std::vector<ExpData*> &edatas,
                                                              const Model &model,
                                                              const bool failfast,
#if defined(_OPENMP)
                                                              int num_threads
#else
                                                              int /* num_threads */
#endif
)
{
    std::vector<std::unique_ptr<ReturnData> > results(edatas.size());
    // is set to true if one simulation fails and we should skip the rest.
    // shared across threads.
    bool skipThrough = false;

#if defined(_OPENMP)
    #pragma omp parallel for num_threads(num_threads)
#endif
    for(int i = 0; i < (int)edatas.size(); ++i) {
        auto mySolver = std::unique_ptr<Solver>(solver.clone());
        auto myModel = std::unique_ptr<Model>(model.clone());

        /* if we fail we need to write empty return datas for the python
         interface */
        if (skipThrough) {
            ConditionContext conditionContext(myModel.get(), edatas[i]);
            results[i] =
                std::unique_ptr<ReturnData>(new ReturnData(solver, model));
        } else {
            results[i] = runAmiciSimulation(*mySolver, edatas[i], *myModel);
        }

        skipThrough |= failfast && results[i]->status < 0;
    }

    return results;
}

} // namespace amici
