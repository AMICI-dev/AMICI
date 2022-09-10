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
#include <map>
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

std::map<int, std::string> simulation_status_to_str_map = {
    {AMICI_RECOVERABLE_ERROR, "AMICI_RECOVERABLE_ERROR"},
    {AMICI_UNRECOVERABLE_ERROR, "AMICI_UNRECOVERABLE_ERROR"},
    {AMICI_TOO_MUCH_WORK, "AMICI_TOO_MUCH_WORK"},
    {AMICI_TOO_MUCH_ACC, "AMICI_TOO_MUCH_ACC"},
    {AMICI_ERR_FAILURE, "AMICI_ERR_FAILURE"},
    {AMICI_CONV_FAILURE, "AMICI_CONV_FAILURE"},
    {AMICI_RHSFUNC_FAIL, "AMICI_RHSFUNC_FAIL"},
    {AMICI_ILL_INPUT, "AMICI_ILL_INPUT"},
    {AMICI_ERROR, "AMICI_ERROR"},
    {AMICI_NO_STEADY_STATE, "AMICI_NO_STEADY_STATE"},
    {AMICI_DAMPING_FACTOR_ERROR, "AMICI_DAMPING_FACTOR_ERROR"},
    {AMICI_SINGULAR_JACOBIAN, "AMICI_SINGULAR_JACOBIAN"},
    {AMICI_NOT_IMPLEMENTED, "AMICI_NOT_IMPLEMENTED"},
    {AMICI_MAX_TIME_EXCEEDED, "AMICI_MAX_TIME_EXCEEDED"},
    {AMICI_SUCCESS, "AMICI_SUCCESS"},
};

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
    auto start_time_total = clock();
    solver.startTimer();

    /* Applies condition-specific model settings and restores them when going
     * out of scope */
    ConditionContext cc1(&model, edata, FixedParameterContext::simulation);

    std::unique_ptr<ReturnData> rdata = std::make_unique<ReturnData>(solver,
                                                                     model);
    if(edata) {
        rdata->id = edata->id;
    }

    std::unique_ptr<SteadystateProblem> preeq {};
    std::unique_ptr<ForwardProblem> fwd {};
    std::unique_ptr<BackwardProblem> bwd {};
    std::unique_ptr<SteadystateProblem> posteq {};
    // tracks whether backwards integration finished without exceptions
    bool bwd_success = true;

    try {
        if (edata && !edata->fixedParametersPreequilibration.empty()) {
            ConditionContext cc2(
                &model, edata, FixedParameterContext::preequilibration
            );

            preeq = std::make_unique<SteadystateProblem>(solver, model);
            preeq->workSteadyStateProblem(solver, model, -1);
        }


        fwd = std::make_unique<ForwardProblem>(edata, &model, &solver,
                                               preeq.get());
        fwd->workForwardProblem();


        if (fwd->getCurrentTimeIteration() < model.nt()) {
            posteq = std::make_unique<SteadystateProblem>(solver, model);
            posteq->workSteadyStateProblem(solver, model,
                                           fwd->getCurrentTimeIteration());
        }


        if (edata && solver.computingASA()) {
            fwd->getAdjointUpdates(model, *edata);
            if (posteq) {
                posteq->getAdjointUpdates(model, *edata);
                posteq->workSteadyStateBackwardProblem(solver, model,
                                                       bwd.get());
            }

            bwd_success = false;

            bwd = std::make_unique<BackwardProblem>(*fwd, posteq.get());
            bwd->workBackwardProblem();

            bwd_success = true;

            if (preeq) {
                ConditionContext cc2(&model, edata,
                                     FixedParameterContext::preequilibration);
                preeq->workSteadyStateBackwardProblem(solver, model,
                                                      bwd.get());
            }
        }

        rdata->status = AMICI_SUCCESS;

    } catch (amici::IntegrationFailure const& ex) {
        if(ex.error_code == AMICI_RHSFUNC_FAIL && solver.timeExceeded()) {
            rdata->status = AMICI_MAX_TIME_EXCEEDED;
            if(rethrow)
                throw;
            warningF("AMICI:simulation",
                     "AMICI forward simulation failed at t = %f: "
                     "Maximum time exceeded.\n",
                     ex.time);
        } else {
            rdata->status = ex.error_code;
            if (rethrow)
                throw;
            warningF("AMICI:simulation",
                     "AMICI forward simulation failed at t = %f:\n%s\n",
                     ex.time,
                     ex.what());

        }
    } catch (amici::IntegrationFailureB const& ex) {
        if(ex.error_code == AMICI_RHSFUNC_FAIL && solver.timeExceeded()) {
            rdata->status = AMICI_MAX_TIME_EXCEEDED;
            if (rethrow)
                throw;
            warningF(
                "AMICI:simulation",
                "AMICI backward simulation failed when trying to solve until "
                "t = %f: Maximum time exceeded.\n",
                ex.time);

        } else {
            rdata->status = ex.error_code;
            if (rethrow)
                throw;
            warningF(
                "AMICI:simulation",
                "AMICI backward simulation failed when trying to solve until t = %f"
                " (see message above):\n%s\n",
                ex.time,
                ex.what());
        }
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

    rdata->cpu_time_total = static_cast<double>(clock() - start_time_total)
                            * 1000.0 / CLOCKS_PER_SEC;

    // verify that reported CPU times are plausible
    gsl_EnsuresDebug(rdata->cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->cpu_timeB <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->preeq_cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->preeq_cpu_timeB <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->posteq_cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->posteq_cpu_timeB <= rdata->cpu_time_total);
    if (!posteq)
        gsl_EnsuresDebug(
            std::is_sorted(rdata->numsteps.begin(), rdata->numsteps.end())
            || rdata->status != AMICI_SUCCESS
        );
    if (!preeq)
        gsl_EnsuresDebug(
            std::is_sorted(rdata->numstepsB.begin(), rdata->numstepsB.end())
            || rdata->status != AMICI_SUCCESS
        );

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

std::string simulation_status_to_str(int status)
{
    try {
        return simulation_status_to_str_map.at(status);
    } catch (std::out_of_range const&) {
        // Missing mapping - terminate if this is a debug build,
        // but show the number if non-debug.
        gsl_ExpectsDebug(false);
        return std::to_string(status);
    }
}

} // namespace amici
