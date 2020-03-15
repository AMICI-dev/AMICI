/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include "amici/amici.h"

#include "amici/backwardproblem.h"
#include "amici/forwardproblem.h"
#include "amici/misc.h"

#include <cvodes/cvodes.h>           //return codes
#include <sundials/sundials_types.h> //realtype

#include <cassert>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
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
static_assert(std::is_same<amici::realtype, realtype>::value,
              "Definition of realtype does not match");

namespace amici {

/** AMICI default application context, kept around for convenience for using
  * amici::runAmiciSimulation or instantiating Solver and Model without special
  * needs.
  */
AmiciApplication defaultContext = AmiciApplication();

std::unique_ptr<ReturnData>
runAmiciSimulation(Solver& solver,
                   const ExpData* edata,
                   Model& model,
                   bool rethrow)
{
    return defaultContext.runAmiciSimulation(solver, edata, model, rethrow);
}

void
printErrMsgIdAndTxt(std::string const& id, std::string const& message)
{
    std::cerr << "[Error] ";
    if (!id.empty()) {
        std::cerr << id << ": ";
    }
    std::cerr << message << std::endl;
}

void
printWarnMsgIdAndTxt(std::string const& id, std::string const& message)
{
    std::cerr << "[Warning] ";
    if (!id.empty()) {
        std::cerr << id << ": ";
    }
    std::cerr << message << std::endl;
}

std::vector<std::unique_ptr<ReturnData>>
runAmiciSimulations(const Solver& solver,
                    const std::vector<ExpData*>& edatas,
                    const Model& model,
                    const bool failfast,
#if defined(_OPENMP)
                    int num_threads
#else
                    int /* num_threads */
#endif
)
{
#if defined(_OPENMP)
    return defaultContext.runAmiciSimulations(
      solver, edatas, model, failfast, num_threads);
#else
    return defaultContext.runAmiciSimulations(solver, edatas, model, failfast, 1);
#endif
}

std::unique_ptr<ReturnData>
AmiciApplication::runAmiciSimulation(Solver& solver,
                                     const ExpData* edata,
                                     Model& model,
                                     bool rethrow)
{
    std::unique_ptr<ReturnData> rdata;

    /* Applies condition-specific model settings and restores them when going
     * out of scope */
    ConditionContext conditionContext(&model, edata);

    try {
        rdata = std::unique_ptr<ReturnData>(new ReturnData(solver, model));

        if (model.nx_solver <= 0) {
            return rdata;
        }

        auto fwd = std::unique_ptr<ForwardProblem>(
          new ForwardProblem(rdata.get(), edata, &model, &solver));
        fwd->workForwardProblem();

        auto bwd =
          std::unique_ptr<BackwardProblem>(new BackwardProblem(fwd.get()));
        bwd->workBackwardProblem();

        rdata->status = AMICI_SUCCESS;
    } catch (amici::IntegrationFailure const& ex) {
        rdata->invalidate(ex.time);
        rdata->status = ex.error_code;
        if (rethrow)
            throw;
        warningF("AMICI:simulation",
                 "AMICI forward simulation failed at t = %f:\n%s\n",
                 ex.time,
                 ex.what());
    } catch (amici::IntegrationFailureB const& ex) {
        rdata->invalidateSLLH();
        rdata->status = ex.error_code;
        if (rethrow)
            throw;
        warningF(
          "AMICI:simulation",
          "AMICI backward simulation failed when trying to solve until t = %f"
          " (see message above):\n%s\n",
          ex.time,
          ex.what());
    } catch (amici::AmiException const& ex) {
        rdata->invalidate(model.t0());
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        warningF("AMICI:simulation",
                 "AMICI simulation failed:\n%s\nError occured in:\n%s",
                 ex.what(),
                 ex.getBacktrace());
    }

    rdata->applyChainRuleFactorToSimulationResults(&model);

    return rdata;
}

std::vector<std::unique_ptr<ReturnData>>
AmiciApplication::runAmiciSimulations(const Solver& solver,
                                      const std::vector<ExpData*>& edatas,
                                      const Model& model,
                                      bool failfast,
#if defined(_OPENMP)
                                      int num_threads
#else
                                      int /* num_threads */
#endif

)
{
    std::vector<std::unique_ptr<ReturnData>> results(edatas.size());
    // is set to true if one simulation fails and we should skip the rest.
    // shared across threads.
    bool skipThrough = false;

#if defined(_OPENMP)
#pragma omp parallel for num_threads(num_threads)
#endif
    for (int i = 0; i < (int)edatas.size(); ++i) {
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

void
AmiciApplication::warningF(const char* identifier, const char* format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    auto str = printfToString(format, argptr);
    va_end(argptr);
    warning(identifier, str);
}

void
AmiciApplication::errorF(const char* identifier, const char* format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    auto str = printfToString(format, argptr);
    va_end(argptr);
    error(identifier, str);
}

int
AmiciApplication::checkFinite(gsl::span<const realtype> array, const char* fun)
{

    for (int idx = 0; idx < (int)array.size(); idx++) {
        if (isNaN(array[idx])) {
            warningF("AMICI:NaN",
                     "AMICI encountered a NaN value at index %i of %i in %s!",
                     idx,
                     (int)array.size(),
                     fun);
            return AMICI_RECOVERABLE_ERROR;
        }
        if (isInf(array[idx])) {
            warningF("AMICI:Inf",
                     "AMICI encountered an Inf value at index %i of %i in %s!",
                     idx,
                     (int)array.size(),
                     fun);
            return AMICI_RECOVERABLE_ERROR;
        }
    }
    return AMICI_SUCCESS;
}

} // namespace amici
