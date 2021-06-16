/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include "amici/amici.h"

#include "amici/steadystateproblem.h"
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
static_assert(amici::AMICI_SUCCESS == CV_SUCCESS,
              "AMICI_SUCCESS != CV_SUCCESS");
static_assert(amici::AMICI_DATA_RETURN == CV_TSTOP_RETURN,
              "AMICI_DATA_RETURN != CV_TSTOP_RETURN");
static_assert(amici::AMICI_ROOT_RETURN == CV_ROOT_RETURN,
              "AMICI_ROOT_RETURN != CV_ROOT_RETURN");
static_assert(amici::AMICI_ILL_INPUT == CV_ILL_INPUT,
              "AMICI_ILL_INPUT != CV_ILL_INPUT");
static_assert(amici::AMICI_NORMAL == CV_NORMAL,
              "AMICI_NORMAL != CV_NORMAL");
static_assert(amici::AMICI_ONE_STEP == CV_ONE_STEP,
              "AMICI_ONE_STEP != CV_ONE_STEP");
static_assert(amici::AMICI_SINGULAR_JACOBIAN == SUNLS_PACKAGE_FAIL_UNREC,
              "AMICI_SINGULAR_JACOBIAN != SUNLS_PACKAGE_FAIL_UNREC");
static_assert(amici::AMICI_SINGULAR_JACOBIAN == SUNLS_PACKAGE_FAIL_UNREC,
              "AMICI_SINGULAR_JACOBIAN != SUNLS_PACKAGE_FAIL_UNREC");
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
    /* Applies condition-specific model settings and restores them when going
     * out of scope */
    ConditionContext cc1(&model, edata, FixedParameterContext::simulation);

    std::unique_ptr<ReturnData> rdata = std::make_unique<ReturnData>(solver,
                                                                     model);

    std::unique_ptr<SteadystateProblem> preeq {};
    std::unique_ptr<ForwardProblem> fwd {};
    std::unique_ptr<BackwardProblem> bwd {};
    std::unique_ptr<SteadystateProblem> posteq {};
    // tracks whether backwards integration finished without exceptions
    bool bwd_success = true;

    try {
        if (solver.getPreequilibration() ||
            (edata && !edata->fixedParametersPreequilibration.empty())) {
            ConditionContext cc2(
                &model, edata, FixedParameterContext::preequilibration
            );

            preeq = std::make_unique<SteadystateProblem>(solver, model);
            preeq->workSteadyStateProblem(&solver, &model, -1);
        }


        fwd = std::make_unique<ForwardProblem>(edata, &model, &solver,
                                               preeq.get());
        fwd->workForwardProblem();


        if (fwd->getCurrentTimeIteration() < model.nt()) {
            posteq = std::make_unique<SteadystateProblem>(solver, model);
            posteq->workSteadyStateProblem(&solver, &model,
                                           fwd->getCurrentTimeIteration());
        }


        if (edata && solver.computingASA()) {
            fwd->getAdjointUpdates(model, *edata);
            if (posteq) {
                posteq->getAdjointUpdates(model, *edata);
                posteq->workSteadyStateBackwardProblem(&solver, &model,
                                                       bwd.get());
            }

            bwd_success = false;

            bwd = std::make_unique<BackwardProblem>(*fwd, posteq.get());
            bwd->workBackwardProblem();

            bwd_success = true;

            if (preeq) {
                ConditionContext cc2(&model, edata,
                                     FixedParameterContext::preequilibration);
                preeq->workSteadyStateBackwardProblem(&solver, &model,
                                                      bwd.get());
            }
        }

        rdata->status = AMICI_SUCCESS;

    } catch (amici::IntegrationFailure const& ex) {
        rdata->status = ex.error_code;
        if (rethrow)
            throw;
        warningF("AMICI:simulation",
                 "AMICI forward simulation failed at t = %f:\n%s\n",
                 ex.time,
                 ex.what());
    } catch (amici::IntegrationFailureB const& ex) {
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
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        warningF("AMICI:simulation",
                 "AMICI simulation failed:\n%s\nError occurred in:\n%s",
                 ex.what(),
                 ex.getBacktrace());
    } catch (std::exception const& ex) {
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        warningF("AMICI:simulation",
                 "AMICI simulation failed:\n%s\n",
                 ex.what());
    }

    rdata->processSimulationObjects(
        preeq.get(), fwd.get(),
        bwd_success ? bwd.get() : nullptr,
        posteq.get(), model, solver, edata);
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
AmiciApplication::warningF(const char* identifier, const char* format, ...) const
{
    va_list argptr;
    va_start(argptr, format);
    auto str = printfToString(format, argptr);
    va_end(argptr);
    warning(identifier, str);
}

void
AmiciApplication::errorF(const char* identifier, const char* format, ...) const
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
                     "AMICI encountered a NaN value at index %i/%i in %s!",
                     idx,
                     (int)array.size()-1,
                     fun);
            return AMICI_RECOVERABLE_ERROR;
        }
        if (isInf(array[idx])) {
            warningF("AMICI:Inf",
                     "AMICI encountered an Inf value at index %i/%i in %s!",
                     idx,
                     (int)array.size()-1,
                     fun);
            return AMICI_RECOVERABLE_ERROR;
        }
    }
    return AMICI_SUCCESS;
}

} // namespace amici
