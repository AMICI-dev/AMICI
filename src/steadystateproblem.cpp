#include "amici/steadystateproblem.h"
#include "amici/backwardproblem.h"
#include "amici/defines.h"
#include "amici/edata.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/newton_solver.h"
#include "amici/solver.h"

#include <cmath>
#include <cvodes/cvodes.h>
#include <memory>
#include <sundials/sundials_dense.h>

namespace amici {

constexpr realtype conv_thresh = 1.0;

/**
 * @brief Assemble the error message to be thrown according to steady state
 * computation status.
 * @param errorString The error string to append to.
 * @param status Status of the steady state computation.
 */
static void
writeErrorString(std::string& errorString, SteadyStateStatus status) {
    switch (status) {
    case SteadyStateStatus::failed_too_long_simulation:
        errorString.append(": System could not be equilibrated.");
        break;
    case SteadyStateStatus::failed_damping:
        errorString.append(": Damping factor reached lower bound.");
        break;
    case SteadyStateStatus::failed_factorization:
        errorString.append(": RHS could not be factorized.");
        break;
    case SteadyStateStatus::failed_convergence:
        errorString.append(": No convergence was achieved.");
        break;
    default:
        errorString.append(".");
        break;
    }
}

/**
 * @brief Compute the weighted root-mean-square norm of xdot.
 *
 * The weights are computed according to x:
 * w_i = 1 / ( rtol * x_i + atol )
 * @param x current state (sx[ip] for sensitivities)
 * @param xdot current rhs (sxdot[ip] for sensitivities)
 * @param mask mask for state variables to include in WRMS norm.
 * Positive value: include; non-positive value: exclude; empty: include all.
 * @param atol absolute tolerance
 * @param rtol relative tolerance
 * @param ewt error weight vector
 * @return root-mean-square norm
 */
realtype getWrmsNorm(
    AmiVector const& x, AmiVector const& xdot, AmiVector const& mask,
    realtype atol, realtype rtol, AmiVector& ewt
) {
    // ewt = x
    N_VAbs(const_cast<N_Vector>(x.getNVector()), ewt.getNVector());
    // ewt *= rtol
    N_VScale(rtol, ewt.getNVector(), ewt.getNVector());
    // ewt += atol
    N_VAddConst(ewt.getNVector(), atol, ewt.getNVector());
    // ewt = 1/ewt (ewt = 1/(rtol*x+atol))
    N_VInv(ewt.getNVector(), ewt.getNVector());

    // wrms = sqrt(sum((xdot/ewt)**2)/n) where n = size of state vector
    if (mask.getLength()) {
        return N_VWrmsNormMask(
            const_cast<N_Vector>(xdot.getNVector()), ewt.getNVector(),
            const_cast<N_Vector>(mask.getNVector())
        );
    }
    return N_VWrmsNorm(
        const_cast<N_Vector>(xdot.getNVector()), ewt.getNVector()
    );
}

realtype WRMSComputer::wrms(AmiVector const& x, AmiVector const& x_ref) {
    return getWrmsNorm(x_ref, x, mask_, atol_, rtol_, ewt_);
}

/**
 * @brief Compute the backward quadratures, which contribute to the
 * gradient (xQB) from the quadrature over the backward state itself (xQ)
 * @param model Model instance.
 * @param yQ vector to be multiplied with dxdotdp
 * @param yQB resulting vector after multiplication
 * @param state Simulation state for which to compute dxdot/dp
 */
void computeQBfromQ(
    Model& model, AmiVector const& yQ, AmiVector& yQB,
    SimulationState const& state
) {
    // Compute the quadrature as the inner product: `yQB = dxdotdp * yQ`.

    // set to zero first, as multiplication adds to existing value
    yQB.zero();
    // yQB += dxdotdp * yQ
    if (model.pythonGenerated) {
        // fill dxdotdp with current values
        auto const& plist = model.getParameterList();
        model.fdxdotdp(state.t, state.x, state.dx);
        model.get_dxdotdp_full().multiply(
            yQB.getNVector(), yQ.getNVector(), plist, true
        );
    } else {
        for (int ip = 0; ip < model.nplist(); ++ip)
            yQB[ip] = dotProd(yQ, model.get_dxdotdp()[ip]);
    }
}

SteadystateProblem::SteadystateProblem(Solver const& solver, Model& model)
    : wrms_computer_x_(
          model.nx_solver, solver.getSunContext(),
          solver.getAbsoluteToleranceSteadyState(),
          solver.getRelativeToleranceSteadyState(),
          AmiVector(model.get_steadystate_mask(), solver.getSunContext())
      )
    , wrms_computer_xQB_(
          model.nplist(), solver.getSunContext(),
          solver.getAbsoluteToleranceQuadratures(),
          solver.getRelativeToleranceQuadratures(), AmiVector()
      )
    , wrms_computer_sx_(
          model.nx_solver, solver.getSunContext(),
          solver.getAbsoluteToleranceSteadyStateSensi(),
          solver.getRelativeToleranceSteadyStateSensi(),
          AmiVector(model.get_steadystate_mask(), solver.getSunContext())
      )
    , xdot_(model.nx_solver, solver.getSunContext())
    , sdx_(model.nx_solver, model.nplist(), solver.getSunContext())
    , xB_(model.nJ * model.nx_solver, solver.getSunContext())
    , xQ_(model.nJ * model.nx_solver, solver.getSunContext())
    , xQB_(model.nplist(), solver.getSunContext())
    , xQBdot_(model.nplist(), solver.getSunContext())
    , state_(
          {.t = INFINITY,
           .x = AmiVector(model.nx_solver, solver.getSunContext()),
           .dx = AmiVector(model.nx_solver, solver.getSunContext()),
           .sx = AmiVectorArray(
               model.nx_solver, model.nplist(), solver.getSunContext()
           ),
           .state = model.getModelState()}
      )
    , newton_solver_(
          NewtonSolver(model, solver.getLinearSolver(), solver.getSunContext())
      )
    , newtons_method_(
          &model, solver.getSunContext(), &newton_solver_,
          solver.getNewtonDampingFactorMode(),
          solver.getNewtonDampingFactorLowerBound(), solver.getNewtonMaxSteps(),
          solver.getNewtonStepSteadyStateCheck()
      )
    , newton_step_conv_(solver.getNewtonStepSteadyStateCheck())
    , check_sensi_conv_(solver.getSensiSteadyStateCheck()) {
    // Check for compatibility of options
    if (solver.getSensitivityMethod() == SensitivityMethod::forward
        && solver.getSensitivityMethodPreequilibration()
               == SensitivityMethod::adjoint
        && solver.getSensitivityOrder() > SensitivityOrder::none)
        throw AmiException(
            "Preequilibration using adjoint sensitivities "
            "is not compatible with using forward "
            "sensitivities during simulation"
        );
    if (solver.getSensitivityMethod() == SensitivityMethod::forward
        && model.getSteadyStateComputationMode()
               == SteadyStateComputationMode::newtonOnly
        && model.getSteadyStateSensitivityMode()
               == SteadyStateSensitivityMode::integrationOnly)
        throw AmiException(
            "For forward sensitivity analysis steady-state "
            "computation mode 'newtonOnly' and steady-state "
            "sensitivity mode 'integrationOnly' are not "
            "compatible as numerical integration of the model "
            "ODEs and corresponding forward sensitivities ODEs "
            "is coupled"
        );
}

void SteadystateProblem::workSteadyStateProblem(
    Solver const& solver, Model& model, int it
) {
    if (model.ne > 0) {
        solver.logger->log(
            LogSeverity::warning, "STEADY_STATE_SIMULATION",
            "Steady-state simulation with events is not supported. "
            "Events will be ignored during pre- and post-equilibration. "
            "This is subject to change."
        );
    }

    initializeForwardProblem(it, solver, model);

    // Compute steady state, track computation time
    CpuTimer cpu_timer;
    findSteadyState(solver, model, it);

    // Check whether state sensitivities still need to be computed.
    if (requires_state_sensitivities(
            model, solver, it, SteadyStateContext::newtonSensi
        )) {
        try {
            // This might still fail, if the Jacobian is singular and
            // simulation did not find a steady state.
            newton_solver_.computeNewtonSensis(state_.sx, model, state_);
        } catch (NewtonFailure const&) {
            throw AmiException(
                "Steady state sensitivity computation failed due "
                "to unsuccessful factorization of RHS Jacobian"
            );
        }
    }
    cpu_time_ = cpu_timer.elapsed_milliseconds();
}

void SteadystateProblem::workSteadyStateBackwardProblem(
    Solver const& solver, Model& model, AmiVector const& xB0, bool is_preeq
) {
    // note that state_ is still set from forward run
    if (is_preeq) {
        if (solver.getSensitivityMethodPreequilibration()
            != SensitivityMethod::adjoint) {
            // if not adjoint mode, there's nothing to do
            return;
        }

        // If we need to reinitialize solver states, this won't work yet.
        if (model.nx_reinit() > 0)
            throw NewtonFailure(
                AMICI_NOT_IMPLEMENTED,
                "Adjoint preequilibration with reinitialization of "
                "non-constant states is not yet implemented. Stopping."
            );

        // only preequilibrations needs a reInit,
        // postequilibration does not
        solver.reInit(state_.t, state_.x, state_.dx);
        solver.updateAndReinitStatesAndSensitivities(&model);
    }

    newton_solver_.reinitialize();
    xB_.copy(xB0);

    // initialize quadratures
    xQ_.zero();
    xQB_.zero();
    xQBdot_.zero();

    // Compute quadratures, track computation time
    CpuTimer cpu_timer;
    computeSteadyStateQuadrature(solver, model);
    cpu_timeB_ = cpu_timer.elapsed_milliseconds();
}

void SteadystateProblem::findSteadyState(
    Solver const& solver, Model& model, int it
) {
    steady_state_status_.resize(3, SteadyStateStatus::not_run);
    // Turn off Newton's method if 'integrationOnly' approach is chosen for
    // steady-state computation or newton_maxsteps is set to 0 or
    // if 'integrationOnly' approach is chosen for sensitivity computation
    // in combination with forward sensitivities approach. The latter is
    // necessary as numerical integration of the model ODEs and corresponding
    // forward sensitivities ODEs is coupled. If 'integrationOnly' approach is
    // chosen for sensitivity computation it is enforced that steady state is
    // computed only by numerical integration as well.
    bool turnOffNewton
        = model.getSteadyStateComputationMode()
              == SteadyStateComputationMode::integrationOnly
          || solver.getNewtonMaxSteps() == 0
          || (solver.getSensitivityOrder() >= SensitivityOrder::first
              && model.getSteadyStateSensitivityMode()
                     == SteadyStateSensitivityMode::integrationOnly
              && ((it == -1
                   && solver.getSensitivityMethodPreequilibration()
                          == SensitivityMethod::forward)
                  || solver.getSensitivityMethod() == SensitivityMethod::forward
              ));

    bool turnOffSimulation = model.getSteadyStateComputationMode()
                             == SteadyStateComputationMode::newtonOnly;

    // First, try to run the Newton solver.
    if (!turnOffNewton)
        findSteadyStateByNewtonsMethod(model, false);

    // Newton solver didn't work, so try to simulate to steady state.
    if (!turnOffSimulation && !checkSteadyStateSuccess())
        steady_state_status_[1]
            = findSteadyStateBySimulation(solver, model, it);

    /* Simulation didn't work, retry the Newton solver from last sim state. */
    if (!turnOffNewton && !turnOffSimulation && !checkSteadyStateSuccess())
        findSteadyStateByNewtonsMethod(model, true);

    // Nothing worked, throw an as informative error as possible.
    if (!checkSteadyStateSuccess())
        handleSteadyStateFailure(
            !turnOffNewton, !turnOffSimulation,
            !turnOffNewton && !turnOffSimulation
        );
}

void SteadystateProblem::findSteadyStateByNewtonsMethod(
    Model& model, bool newton_retry
) {
    int stage = newton_retry ? 2 : 0;
    try {
        updateRightHandSide(model);
        newtons_method_.run(xdot_, state_, wrms_computer_x_);
        steady_state_status_[stage] = SteadyStateStatus::success;
    } catch (NewtonFailure const& ex) {
        switch (ex.error_code) {
        case AMICI_TOO_MUCH_WORK:
            steady_state_status_[stage] = SteadyStateStatus::failed_convergence;
            break;
        case AMICI_NO_STEADY_STATE:
            steady_state_status_[stage]
                = SteadyStateStatus::failed_too_long_simulation;
            break;
        case AMICI_SINGULAR_JACOBIAN:
            steady_state_status_[stage]
                = SteadyStateStatus::failed_factorization;
            break;
        case AMICI_DAMPING_FACTOR_ERROR:
            steady_state_status_[stage] = SteadyStateStatus::failed_damping;
            break;
        default:
            steady_state_status_[stage] = SteadyStateStatus::failed;
            break;
        }
    }
    numsteps_.at(stage) = newtons_method_.get_num_steps();
    wrms_ = newtons_method_.get_wrms();
    flagUpdatedState();
}

SteadyStateStatus SteadystateProblem::findSteadyStateBySimulation(
    Solver const& solver, Model& model, int it
) {
    try {
        if (it < 0) {
            // Preequilibration? -> Create a new solver instance for simulation
            bool integrateSensis = requires_state_sensitivities(
                model, solver, it, SteadyStateContext::solverCreation
            );
            auto newtonSimSolver = createSteadystateSimSolver(
                solver, model, integrateSensis, false
            );
            runSteadystateSimulationFwd(*newtonSimSolver, model);
        } else {
            // Solver was already created, use this one
            runSteadystateSimulationFwd(solver, model);
        }
        return SteadyStateStatus::success;
    } catch (IntegrationFailure const& ex) {
        switch (ex.error_code) {
        case AMICI_TOO_MUCH_WORK:
            if (model.logger)
                model.logger->log(
                    LogSeverity::debug, "EQUILIBRATION_FAILURE",
                    "AMICI equilibration exceeded maximum number of"
                    " integration steps at t=%g.",
                    ex.time
                );
            return SteadyStateStatus::failed_convergence;
        case AMICI_RHSFUNC_FAIL:
            if (model.logger)
                model.logger->log(
                    LogSeverity::debug, "EQUILIBRATION_FAILURE",
                    "AMICI equilibration was stopped after exceedingly"
                    " long simulation time at t=%g.",
                    ex.time
                );
            return SteadyStateStatus::failed_too_long_simulation;
        default:
            if (model.logger)
                model.logger->log(
                    LogSeverity::debug, "OTHER",
                    "AMICI equilibration failed at t=%g.", ex.time
                );
            return SteadyStateStatus::failed;
        }
    } catch (AmiException const& ex) {
        if (model.logger)
            model.logger->log(
                LogSeverity::debug, "OTHER", "AMICI equilibration failed: %s",
                ex.what()
            );
        return SteadyStateStatus::failed;
    }
}

void SteadystateProblem::initializeForwardProblem(
    int it, Solver const& solver, Model& model
) {
    newton_solver_.reinitialize();
    // Process solver handling for pre- or postequilibration.
    if (it == -1) {
        // The solver was not run before, set up everything.
        auto roots_found = std::vector<int>(model.ne, 0);
        model.initialize(
            state_.x, state_.dx, state_.sx, sdx_,
            solver.getSensitivityOrder() >= SensitivityOrder::first, roots_found
        );
        state_.t = model.t0();
        solver.setup(state_.t, &model, state_.x, state_.dx, state_.sx, sdx_);
    } else {
        // The solver was run before, extract current state from solver.
        solver.writeSolution(&state_.t, state_.x, state_.dx, state_.sx, xQ_);
    }

    // overwrite starting timepoint
    if (it < 1) {
        // No previous time point computed, set t = t0
        state_.t = model.t0();
    } else {
        // Carry on simulating from last point
        state_.t = model.getTimepoint(it - 1);
    }

    state_.state = model.getModelState();
    flagUpdatedState();
}

void SteadystateProblem::computeSteadyStateQuadrature(
    Solver const& solver, Model& model
) {
    // This routine computes the quadratures:
    //     xQB = Integral[ xB(x(t), t, p) * dxdot/dp(x(t), t, p) | dt ]
    // As we're in steady state, we have x(t) = x_ss (x_steadystate), hence
    //     xQB = Integral[ xB(x_ss, t, p) | dt ] * dxdot/dp(x_ss, t, p)
    // We therefore compute the integral over xB first and then do a
    // matrix-vector multiplication.

    auto sensitivityMode = model.getSteadyStateSensitivityMode();

    // Try to compute the analytical solution for quadrature algebraically
    if (sensitivityMode == SteadyStateSensitivityMode::newtonOnly
        || sensitivityMode
               == SteadyStateSensitivityMode::integrateIfNewtonFails)
        getQuadratureByLinSolve(model);

    // Perform simulation if necessary
    if (sensitivityMode == SteadyStateSensitivityMode::integrationOnly
        || (sensitivityMode
                == SteadyStateSensitivityMode::integrateIfNewtonFails
            && !hasQuadrature()))
        getQuadratureBySimulation(solver, model);

    // If the analytic solution and integration did not work, throw
    if (!hasQuadrature())
        throw AmiException(
            "Steady state backward computation failed: Linear "
            "system could not be solved (possibly due to singular Jacobian), "
            "and numerical integration did not equilibrate within maxsteps"
        );
}

void SteadystateProblem::getQuadratureByLinSolve(Model& model) {
    // Computes the integral over the adjoint state xB:
    // If the Jacobian has full rank, this has an analytical solution, since
    //   d/dt[ xB(t) ] = JB^T(x(t), p) xB(t) = JB^T(x_ss, p) xB(t)
    // This linear ODE system with time-constant matrix has the solution
    //   xB(t) = exp( t * JB^T(x_ss, p) ) * xB(0)
    // This integral xQ over xB is given as the solution of
    //   JB^T(x_ss, p) * xQ = xB(0)
    // So we first try to solve the linear system, if possible.

    // copy content of xB into vector with integral
    xQ_.copy(xB_);

    // try to solve the linear system
    try {
        // compute integral over xB and write to xQ
        newton_solver_.prepareLinearSystemB(model, state_);
        newton_solver_.solveLinearSystem(xQ_);
        // Compute the quadrature as the inner product xQ * dxdotdp
        computeQBfromQ(model, xQ_, xQB_, state_);
        hasQuadrature_ = true;

        // Finalize by setting adjoint state to zero (its steady state)
        xB_.zero();
    } catch (NewtonFailure const&) {
        hasQuadrature_ = false;
    }
}

void SteadystateProblem::getQuadratureBySimulation(
    Solver const& solver, Model& model
) {
    // If the Jacobian is singular, the integral over xB must be computed
    // by usual integration over time, but simplifications can be applied:
    // x is not time-dependent, no forward trajectory is needed.

    // Set starting timepoint for the simulation solver
    state_.t = model.t0();
    // xQ was written in getQuadratureByLinSolve() -> set to zero
    xQ_.zero();

    auto simSolver = createSteadystateSimSolver(solver, model, false, true);

    // perform integration and quadrature
    try {
        runSteadystateSimulationBwd(*simSolver, model);
        hasQuadrature_ = true;
    } catch (NewtonFailure const&) {
        hasQuadrature_ = false;
    }
}

[[noreturn]] void SteadystateProblem::handleSteadyStateFailure(
    bool tried_newton_1, bool tried_simulation, bool tried_newton_2
) const {
    // Throw error message according to error codes
    std::string errorString = "Steady state computation failed.";
    if (tried_newton_1) {
        errorString.append(" First run of Newton solver failed");
        writeErrorString(errorString, steady_state_status_[0]);
    }
    if (tried_simulation) {
        errorString.append(" Simulation to steady state failed");
        writeErrorString(errorString, steady_state_status_[1]);
    }
    if (tried_newton_2) {
        errorString.append(" Second run of Newton solver failed");
        writeErrorString(errorString, steady_state_status_[2]);
    }
    throw AmiException(errorString.c_str());
}

bool SteadystateProblem::requires_state_sensitivities(
    Model const& model, Solver const& solver, int it, SteadyStateContext context
) const {
    // We need to check whether we need to compute forward sensitivities.
    // Depending on the situation (pre-/postequilibration) and the solver
    // settings, the logic may be involved and is handled here.
    // Most of these boolean operation could be simplified. However,
    // clarity is more important than brevity.

    // Are we running in preequilibration (and hence create)?
    bool preequilibration = (it == -1);

    // Did we already compute forward sensitivities?
    bool forwardSensisAlreadyComputed
        = solver.getSensitivityOrder() >= SensitivityOrder::first
          && steady_state_status_[1] == SteadyStateStatus::success
          && (model.getSteadyStateSensitivityMode()
                  == SteadyStateSensitivityMode::integrationOnly
              || model.getSteadyStateSensitivityMode()
                     == SteadyStateSensitivityMode::integrateIfNewtonFails);

    bool simulationStartedInSteadystate
        = steady_state_status_[0] == SteadyStateStatus::success
          && numsteps_[0] == 0;

    // Do we need forward sensis for postequilibration?
    bool needForwardSensisPosteq
        = !preequilibration && !forwardSensisAlreadyComputed
          && solver.getSensitivityOrder() >= SensitivityOrder::first
          && solver.getSensitivityMethod() == SensitivityMethod::forward;

    // Do we need forward sensis for preequilibration?
    bool needForwardSensisPreeq
        = preequilibration && !forwardSensisAlreadyComputed
          && solver.getSensitivityMethodPreequilibration()
                 == SensitivityMethod::forward
          && solver.getSensitivityOrder() >= SensitivityOrder::first;

    // Do we need to do the linear system solve to get forward sensitivities?
    bool needForwardSensisNewton
        = (needForwardSensisPreeq || needForwardSensisPosteq)
          && !simulationStartedInSteadystate;

    // When we're creating a new solver object
    bool needForwardSensiAtCreation
        = needForwardSensisPreeq
          && (model.getSteadyStateSensitivityMode()
                  == SteadyStateSensitivityMode::integrationOnly
              || model.getSteadyStateSensitivityMode()
                     == SteadyStateSensitivityMode::integrateIfNewtonFails);

    // Check if we need to store sensitivities
    switch (context) {
    case SteadyStateContext::newtonSensi:
        return needForwardSensisNewton;

    case SteadyStateContext::sensiStorage:
        return needForwardSensisNewton || forwardSensisAlreadyComputed
               || simulationStartedInSteadystate;

    case SteadyStateContext::solverCreation:
        return needForwardSensiAtCreation;

    default:
        throw AmiException(
            "Requested invalid context in sensitivity "
            "processing during steady state computation"
        );
    }
}

realtype SteadystateProblem::getWrmsState(Model& model) {
    updateRightHandSide(model);

    if (newton_step_conv_) {
        newtons_method_.compute_step(xdot_, state_);
        return wrms_computer_x_.wrms(newtons_method_.get_delta(), state_.x);
    }

    return wrms_computer_x_.wrms(xdot_, state_.x);
}

realtype SteadystateProblem::getWrmsFSA(Model& model) {
    // Forward sensitivities: Compute weighted error norm for their RHS
    realtype wrms = 0.0;

    // we don't need to call prepareLinearSystem in this function, since it was
    // already computed in the preceding getWrms call and both equations have
    // the same Jacobian.

    xdot_updated_ = false;
    for (int ip = 0; ip < model.nplist(); ++ip) {
        model.fsxdot(
            state_.t, state_.x, state_.dx, ip, state_.sx[ip], state_.dx, xdot_
        );
        if (newton_step_conv_) {
            newton_solver_.solveLinearSystem(xdot_);
        }
        wrms = wrms_computer_sx_.wrms(xdot_, state_.sx[ip]);
        // ideally this function would report the maximum of all wrms over
        // all ip, but for practical purposes we can just report the wrms for
        // the first ip where we know that the convergence threshold is not
        // satisfied.
        if (wrms > conv_thresh)
            break;
    }
    // Just report the parameter for the last `ip`. The value doesn't matter -
    // it's only important that all of them satisfy the convergence threshold.
    return wrms;
}

bool SteadystateProblem::checkSteadyStateSuccess() const {
    // Did one of the attempts yield a steady state?
    return std::ranges::any_of(
        steady_state_status_, [](SteadyStateStatus status
                              ) { return status == SteadyStateStatus::success; }
    );
}

void SteadystateProblem::runSteadystateSimulationFwd(
    Solver const& solver, Model& model
) {
    if (model.nx_solver == 0)
        return;

    // Do we also have to check for convergence of sensitivities?
    SensitivityMethod sensitivity_method = SensitivityMethod::none;
    if (solver.getSensitivityOrder() > SensitivityOrder::none
        && solver.getSensitivityMethod() == SensitivityMethod::forward) {
        sensitivity_method = SensitivityMethod::forward;
    }
    // If forward sensitivity computation by simulation is disabled,
    // disable forward sensitivity integration in the solver.
    // Sensitivities will be computed by newtonsolver.computeNewtonSensis then.
    if (model.getSteadyStateSensitivityMode()
        == SteadyStateSensitivityMode::newtonOnly) {
        solver.switchForwardSensisOff();
        sensitivity_method = SensitivityMethod::none;
    }

    int& sim_steps = numsteps_.at(1);
    int convergence_check_frequency = newton_step_conv_ ? 25 : 1;

    while (true) {
        if (sim_steps % convergence_check_frequency == 0) {
            // Check for convergence (already before simulation, since we might
            // start in steady state)
            wrms_ = getWrmsState(model);
            if (wrms_ < conv_thresh) {
                if (check_sensi_conv_
                    && sensitivity_method == SensitivityMethod::forward) {
                    updateSensiSimulation(solver);
                    // getWrms needs to be called before getWrmsFSA
                    // such that the linear system is prepared for newton-type
                    // convergence check
                    if (getWrmsFSA(model) < conv_thresh)
                        break; // converged
                } else {
                    break; // converged
                }
            }
        }

        // check for maxsteps
        if (sim_steps >= solver.getMaxSteps()) {
            throw IntegrationFailure(AMICI_TOO_MUCH_WORK, state_.t);
        }

        // increase counter
        sim_steps++;

        // One step of ODE integration
        // Reason for tout specification:
        // * max with 1 ensures the correct direction
        //  (any positive value would do)
        // * multiplication with 10 ensures nonzero difference and should
        //   ensure stable computation.
        // The value is not important for AMICI_ONE_STEP mode, only the
        // direction w.r.t. current t.
        solver.step(std::max(state_.t, 1.0) * 10);

        solver.writeSolution(&state_.t, state_.x, state_.dx, state_.sx, xQ_);
        flagUpdatedState();
    }

    // if check_sensi_conv_ is deactivated, we still have to update sensis
    if (sensitivity_method == SensitivityMethod::forward)
        updateSensiSimulation(solver);
}

void SteadystateProblem::runSteadystateSimulationBwd(
    Solver const& solver, Model& model
) {
    if (model.nx_solver == 0)
        return;

    if (newton_step_conv_) {
        throw NewtonFailure(
            AMICI_NOT_IMPLEMENTED,
            "Newton type convergence check is not implemented for adjoint "
            "steady state computations. Stopping."
        );
    }

    int& sim_steps = numstepsB_;

    int convergence_check_frequency = newton_step_conv_ ? 25 : 1;
    auto max_steps = (solver.getMaxStepsBackwardProblem() > 0)
                         ? solver.getMaxStepsBackwardProblem()
                         : solver.getMaxSteps() * 100;

    while (true) {
        if (sim_steps % convergence_check_frequency == 0) {
            // Check for convergence (already before simulation, since we might
            // start in steady state)

            // In the adjoint case, only xQB contributes to the gradient, the
            // exact steadystate is less important, as xB = xQdot may even not
            // converge to zero at all. So we need xQBdot, hence compute xQB
            // first.
            computeQBfromQ(model, xQ_, xQB_, state_);
            computeQBfromQ(model, xB_, xQBdot_, state_);
            wrms_ = wrms_computer_xQB_.wrms(xQBdot_, xQB_);
            if (wrms_ < conv_thresh) {
                break; // converged
            }
        }

        // check for maxsteps
        if (sim_steps >= max_steps) {
            throw IntegrationFailureB(AMICI_TOO_MUCH_WORK, state_.t);
        }

        // increase counter
        sim_steps++;

        // One step of ODE integration
        // Reason for tout specification:
        // * max with 1 ensures the correct direction
        //  (any positive value would do)
        // * multiplication with 10 ensures nonzero difference and should
        //   ensure stable computation.
        // The value is not important for AMICI_ONE_STEP mode, only the
        // direction w.r.t. current t.
        solver.step(std::max(state_.t, 1.0) * 10);

        solver.writeSolution(&state_.t, xB_, state_.dx, state_.sx, xQ_);
    }
}

std::unique_ptr<Solver> SteadystateProblem::createSteadystateSimSolver(
    Solver const& solver, Model& model, bool forwardSensis, bool backward
) const {
    switch (solver.getLinearSolver()) {
    case LinearSolver::dense:
    case LinearSolver::KLU:
        break;
    default:
        throw AmiException("Invalid solver for steady state simulation");
    }

    auto sim_solver = std::unique_ptr<Solver>(solver.clone());

    sim_solver->logger = solver.logger;

    // do we need sensitivities?
    if (forwardSensis) {
        // need forward to compute sx0
        sim_solver->setSensitivityMethod(SensitivityMethod::forward);
    } else {
        sim_solver->setSensitivityMethod(SensitivityMethod::none);
        sim_solver->setSensitivityOrder(SensitivityOrder::none);
    }
    // use x and sx as dummies for dx and sdx
    // (they won't get touched in a CVodeSolver)
    if (backward) {
        sim_solver->setup(model.t0(), &model, xB_, xB_, state_.sx, sdx_);
        sim_solver->setupSteadystate(
            model.t0(), &model, state_.x, state_.dx, xB_, xB_, xQ_
        );
    } else {
        sim_solver->setup(
            model.t0(), &model, state_.x, state_.dx, state_.sx, sdx_
        );
    }

    return sim_solver;
}

void SteadystateProblem::getAdjointUpdates(
    Model& model, ExpData const& edata, std::vector<realtype>& dJydx
) {
    for (int it = 0; it < model.nt(); it++) {
        if (std::isinf(model.getTimepoint(it))) {
            model.getAdjointStateObservableUpdate(
                slice(dJydx, it, model.nx_solver * model.nJ), it, state_.x,
                edata
            );
        }
    }
}

void SteadystateProblem::flagUpdatedState() {
    xdot_updated_ = false;
    sensis_updated_ = false;
}

void SteadystateProblem::updateSensiSimulation(Solver const& solver) {
    if (sensis_updated_)
        return;
    state_.sx = solver.getStateSensitivity(state_.t);
    sensis_updated_ = true;
}

void SteadystateProblem::updateRightHandSide(Model& model) {
    if (xdot_updated_)
        return;
    model.fxdot(state_.t, state_.x, state_.dx, xdot_);
    xdot_updated_ = true;
}

NewtonsMethod::NewtonsMethod(
    gsl::not_null<Model*> model, SUNContext sunctx,
    gsl::not_null<NewtonSolver*> solver,
    NewtonDampingFactorMode damping_factor_mode,
    realtype damping_factor_lower_bound, int max_steps, bool check_delta
)
    : model_(model)
    , max_steps_(max_steps)
    , damping_factor_mode_(damping_factor_mode)
    , damping_factor_lower_bound_(damping_factor_lower_bound)
    , check_delta_(check_delta)
    , solver_(solver)
    , delta_(model->nx_solver, sunctx)
    , delta_old_(model->nx_solver, sunctx)
    , x_old_(model->nx_solver, sunctx) {}

void NewtonsMethod::run(
    AmiVector& xdot, SimulationState& state, WRMSComputer& wrms_computer
) {
    i_step = 0;

    if (model_->nx_solver == 0) {
        wrms_ = 0.0;
        return;
    }

    wrms_ = INFINITY;
    delta_.zero();

    // The Newton step size.
    double gamma{1.0};
    bool update_direction = true;

    wrms_ = compute_wrms(xdot, state, wrms_computer);
    bool converged = has_converged(xdot, state, wrms_computer);

    // Whether the step was successful
    bool step_successful = true;

    while (!converged && i_step < max_steps_) {
        if (step_successful) {
            // If new residuals are smaller than the old ones, update state
            x_old_.copy(state.x);
        }

        // If Newton steps are necessary, compute the initial search
        // direction
        if (update_direction) {
            // compute the next step if not already done during the previous
            // delta-convergence check
            if (!check_delta_) {
                compute_step(xdot, state);
            };

            // we store delta_ here as later convergence checks may update
            // it
            delta_old_.copy(delta_);
        }

        // Try step with new gamma_/delta_, evaluate rhs
        // x = x_old + delta_[old_] * gamma
        linearSum(
            1.0, x_old_, gamma, update_direction ? delta_ : delta_old_, state.x
        );
        model_->fxdot(state.t, state.x, state.dx, xdot);

        realtype wrms_tmp = compute_wrms(xdot, state, wrms_computer);
        step_successful = wrms_tmp < wrms_;
        if (step_successful) {
            wrms_ = wrms_tmp;
            converged = has_converged(xdot, state, wrms_computer);
        }

        update_direction = update_damping_factor(step_successful, gamma);
        ++i_step;
    }

    if (!converged)
        throw NewtonFailure(AMICI_TOO_MUCH_WORK, "applyNewtonsMethod");
}

void NewtonsMethod::compute_step(
    AmiVector const& xdot, SimulationState const& state
) {
    delta_.copy(xdot);
    solver_->getStep(delta_, *model_, state);
}

bool NewtonsMethod::update_damping_factor(bool step_successful, double& gamma) {
    if (damping_factor_mode_ != NewtonDampingFactorMode::on)
        return true;

    if (step_successful) {
        gamma = fmin(1.0, 2.0 * gamma);
    } else {
        gamma /= 4.0;
    }

    if (gamma < damping_factor_lower_bound_) {
        throw NewtonFailure(
            AMICI_DAMPING_FACTOR_ERROR,
            "Newton solver failed: the damping factor "
            "reached its lower bound"
        );
    }
    return step_successful;
}

realtype NewtonsMethod::compute_wrms(
    AmiVector const& xdot, SimulationState const& state,
    WRMSComputer& wrms_computer
) {
    if (check_delta_) {
        compute_step(xdot, state);
        return wrms_computer.wrms(delta_, state.x);
    } else {
        return wrms_computer.wrms(xdot, state.x);
    }
}

bool NewtonsMethod::has_converged(
    AmiVector& xdot, SimulationState& state, WRMSComputer& wrms_computer
) {
    // pre-check convergence
    if (wrms_ >= conv_thresh)
        return false;

    if (!model_->get_any_state_nonnegative()) {
        // no constraints to check for
        return true;
    }

    // Ensure state positivity if requested,
    // and repeat the convergence check if necessary
    auto nonnegative = model_->getStateIsNonNegative();
    Expects(nonnegative.size() == state.x.getVector().size());
    auto state_modified = false;
    for (int ix = 0; ix < state.x.getLength(); ix++) {
        if (state.x[ix] < 0.0 && nonnegative[ix]) {
            state.x[ix] = 0.0;
            state_modified = true;
        }
    }
    if (!state_modified)
        return true;

    model_->fxdot(state.t, state.x, state.dx, xdot);
    wrms_ = compute_wrms(xdot, state, wrms_computer);

    return wrms_ < conv_thresh;
}

} // namespace amici
