/**
 * @file   amici.cpp
 * @brief  core routines for integration
 */

#include "amici/amici.h"

#include "amici/backwardproblem.h"
#include "amici/forwardproblem.h"
#include "amici/logging.h"
#include "amici/steadystateproblem.h"

#include <sundials/sundials_types.h> //sunrealtype

#include <map>
#include <memory>
#include <type_traits>

static_assert(
    amici::AMICI_SINGULAR_JACOBIAN == SUN_ERR_EXT_FAIL,
    "AMICI_SINGULAR_JACOBIAN != SUN_ERR_EXT_FAIL"
);
static_assert(
    std::is_same_v<amici::realtype, sunrealtype>,
    "Definition of realtype does not match"
);

namespace amici {

std::map<int, std::string> simulation_status_to_str_map = {
    {AMICI_RECOVERABLE_ERROR, "AMICI_RECOVERABLE_ERROR"},
    {AMICI_UNRECOVERABLE_ERROR, "AMICI_UNRECOVERABLE_ERROR"},
    {AMICI_TOO_MUCH_WORK, "AMICI_TOO_MUCH_WORK"},
    {AMICI_TOO_MUCH_ACC, "AMICI_TOO_MUCH_ACC"},
    {AMICI_ERR_FAILURE, "AMICI_ERR_FAILURE"},
    {AMICI_CONV_FAILURE, "AMICI_CONV_FAILURE"},
    {AMICI_FIRST_RHSFUNC_ERR, "AMICI_FIRST_RHSFUNC_ERR"},
    {AMICI_CONSTR_FAIL, "AMICI_CONSTR_FAIL"},
    {AMICI_CVODES_CONSTR_FAIL, "AMICI_CVODES_CONSTR_FAIL"},
    {AMICI_IDAS_CONSTR_FAIL, "AMICI_IDAS_CONSTR_FAIL"},
    {AMICI_RHSFUNC_FAIL, "AMICI_RHSFUNC_FAIL"},
    {AMICI_ILL_INPUT, "AMICI_ILL_INPUT"},
    {AMICI_ERROR, "AMICI_ERROR"},
    {AMICI_NO_STEADY_STATE, "AMICI_NO_STEADY_STATE"},
    {AMICI_DAMPING_FACTOR_ERROR, "AMICI_DAMPING_FACTOR_ERROR"},
    {AMICI_SINGULAR_JACOBIAN, "AMICI_SINGULAR_JACOBIAN"},
    {AMICI_NOT_IMPLEMENTED, "AMICI_NOT_IMPLEMENTED"},
    {AMICI_MAX_TIME_EXCEEDED, "AMICI_MAX_TIME_EXCEEDED"},
    {AMICI_SUCCESS, "AMICI_SUCCESS"},
    {AMICI_NOT_RUN, "AMICI_NOT_RUN"},
    {AMICI_LSETUP_FAIL, "AMICI_LSETUP_FAIL"},
    {AMICI_FIRST_QRHSFUNC_ERR, "AMICI_FIRST_QRHSFUNC_ERR"},
    {AMICI_WARNING, "AMICI_WARNING"},
    {AMICI_BAD_T, "AMICI_BAD_T"},
    {AMICI_BAD_DKY, "AMICI_BAD_DKY"},
    {AMICI_FIRST_SRHSFUNC_ERR, "AMICI_FIRST_SRHSFUNC_ERR"},
    {AMICI_SRHSFUNC_FAIL, "AMICI_SRHSFUNC_FAIL"},
    {AMICI_REPTD_SRHSFUNC_ERR, "AMICI_REPTD_SRHSFUNC_ERR"},
    {AMICI_UNREC_SRHSFUNC_ERR, "AMICI_UNREC_SRHSFUNC_ERR"},
    {AMICI_RTFUNC_FAIL, "AMICI_RTFUNC_FAIL"},
    {AMICI_LINESEARCH_FAIL, "AMICI_LINESEARCH_FAIL"},
};

std::unique_ptr<ReturnData> runAmiciSimulation(
    Solver& solver, ExpData const* edata, Model& model, bool const rethrow
) {
    // create a temporary logger instance for Solver and Model to capture
    // messages from only this simulation
    Logger logger;
    solver.logger = &logger;
    model.logger = &logger;
    // prevent dangling pointer
    auto _ = gsl::finally([&solver, &model] {
        solver.logger = model.logger = nullptr;
    });

    CpuTimer cpu_timer;
    solver.startTimer();

    // Applies condition-specific model settings and restores them when going
    // out of scope. (This also sets `plist`, which is required for initializing
    // ReturnData below.)
    ConditionContext cc1(&model, edata, FixedParameterContext::simulation);

    auto rdata = std::make_unique<ReturnData>(solver, model);
    if (edata) {
        rdata->id = edata->id;
    }

    std::unique_ptr<ForwardProblem> fwd{};
    std::unique_ptr<BackwardProblem> bwd{};
    // tracks whether backwards integration finished without exceptions
    bool bwd_success = true;

    try {
        fwd = std::make_unique<ForwardProblem>(edata, &model, &solver);
        fwd->workForwardProblem();

        if (edata && solver.computingASA()) {
            bwd_success = false; // NOLINT

            bwd = std::make_unique<BackwardProblem>(*fwd);
            bwd->workBackwardProblem();

            bwd_success = true;
        }

        rdata->status = AMICI_SUCCESS;
    } catch (IntegrationFailure const& ex) {
        if (ex.error_code == AMICI_RHSFUNC_FAIL && solver.timeExceeded()) {
            rdata->status = AMICI_MAX_TIME_EXCEEDED;
            if (rethrow)
                throw;
            logger.log(
                LogSeverity::error, "MAXTIME_EXCEEDED",
                "AMICI forward simulation failed at t = %g: "
                "Maximum time exceeded in forward solve.",
                ex.time
            );
        } else {
            rdata->status = ex.error_code;
            if (rethrow)
                throw;
            logger.log(
                LogSeverity::error, "FORWARD_FAILURE",
                "AMICI forward simulation failed at t = %g: %s", ex.time,
                ex.what()
            );
        }
    } catch (IntegrationFailureB const& ex) {
        if (ex.error_code == AMICI_RHSFUNC_FAIL && solver.timeExceeded()) {
            rdata->status = AMICI_MAX_TIME_EXCEEDED;
            if (rethrow)
                throw;
            logger.log(
                LogSeverity::error, "MAXTIME_EXCEEDED",
                "AMICI backward simulation failed when trying to solve until "
                "t = %g: Maximum time exceeded in backward solve.",
                ex.time
            );

        } else {
            rdata->status = ex.error_code;
            if (rethrow)
                throw;
            logger.log(
                LogSeverity::error, "BACKWARD_FAILURE",
                "AMICI backward simulation failed when trying to solve until t "
                "= %g"
                " (check debug logs for details): %s",
                ex.time, ex.what()
            );
        }
    } catch (AmiException const& ex) {
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        logger.log(
            LogSeverity::error, "OTHER", "AMICI simulation failed: %s",
            ex.what()
        );
#ifndef NDEBUG
        logger.log(
            LogSeverity::debug, "BACKTRACE",
            "The previous error occurred at:\n%s", ex.getBacktrace()
        );
#endif
    } catch (std::exception const& ex) {
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        logger.log(
            LogSeverity::error, "OTHER", "AMICI simulation failed: %s",
            ex.what()
        );
    }

    try {
        rdata->processSimulationObjects(
            fwd.get(), bwd_success ? bwd.get() : nullptr, model, solver, edata
        );
    } catch (std::exception const& ex) {
        rdata->status = AMICI_ERROR;
        if (rethrow)
            throw;
        logger.log(
            LogSeverity::error, "OTHER", "AMICI simulation failed: %s",
            ex.what()
        );
    }

    rdata->t_last = solver.gett();
    rdata->cpu_time_total = cpu_timer.elapsed_milliseconds();

    // verify that reported CPU times are plausible
    gsl_EnsuresDebug(rdata->cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->cpu_timeB <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->preeq_cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->preeq_cpu_timeB <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->posteq_cpu_time <= rdata->cpu_time_total);
    gsl_EnsuresDebug(rdata->posteq_cpu_timeB <= rdata->cpu_time_total);
    if (fwd && !fwd->getPostequilibrationProblem())
        gsl_EnsuresDebug(
            std::ranges::is_sorted(rdata->numsteps)
            || rdata->status != AMICI_SUCCESS
        );
    if (fwd && !fwd->getPreequilibrationProblem())
        gsl_EnsuresDebug(
            std::ranges::is_sorted(rdata->numstepsB)
            || rdata->status != AMICI_SUCCESS
        );

    rdata->messages = logger.items;

    return rdata;
}

std::vector<std::unique_ptr<ReturnData>> runAmiciSimulations(
    Solver const& solver, std::vector<ExpData*> const& edatas,
    Model const& model, bool const failfast,
#if defined(_OPENMP)
    int num_threads
#else
    int /* num_threads */
#endif

) {
    std::vector<std::unique_ptr<ReturnData>> results(edatas.size());
    // is set to true if one simulation fails and we should skip the rest.
    // shared across threads.
    bool skipThrough = false;

#if defined(_OPENMP)
#pragma omp parallel for num_threads(num_threads)
#endif
    for (int i = 0; i < (int)edatas.size(); ++i) {
        // must catch exceptions in parallel section to avoid termination
        try {
            auto mySolver = std::unique_ptr<Solver>(solver.clone());
            auto myModel = std::unique_ptr<Model>(model.clone());

            /* if we fail we need to write empty return datas for the python
             interface */
            if (skipThrough) {
                ConditionContext conditionContext(myModel.get(), edatas[i]);
                results[i] = std::make_unique<ReturnData>(solver, model);
            } else {
                results[i] = runAmiciSimulation(*mySolver, edatas[i], *myModel);
            }
        } catch (std::exception const& ex) {
            results[i] = std::make_unique<ReturnData>(solver, model);
            results[i]->status = AMICI_ERROR;
            results[i]->messages.emplace_back(
                LogSeverity::error, "OTHER", ex.what()
            );
        }

        skipThrough |= failfast && results[i]->status < 0;
    }

    return results;
}

std::string simulation_status_to_str(int const status) {
    try {
        return simulation_status_to_str_map.at(status);
    } catch (std::out_of_range const&) {
        // Missing mapping - terminate if this is a debug build,
        // but show the number if non-debug.
        fprintf(stderr, "Unknown simulation status: %d\n", status);
        gsl_ExpectsDebug(false);
        return std::to_string(status);
    }
}

} // namespace amici
