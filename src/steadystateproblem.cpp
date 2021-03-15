#include "amici/steadystateproblem.h"
#include "amici/defines.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/solver_cvodes.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
#include "amici/backwardproblem.h"
#include "amici/newton_solver.h"
#include "amici/misc.h"

#include <cmath>
#include <cstring>
#include <ctime>
#include <sundials/sundials_dense.h>
#include <memory>
#include <cvodes/cvodes.h>

namespace amici {

SteadystateProblem::SteadystateProblem(const Solver &solver, const Model &model)
    : delta_(model.nx_solver), ewt_(model.nx_solver), ewtQB_(model.nplist()),
      rel_x_newton_(model.nx_solver), x_newton_(model.nx_solver),
      x_(model.nx_solver), x_old_(model.nx_solver), dx_(model.nx_solver),
      xdot_(model.nx_solver), xdot_old_(model.nx_solver),
      sx_(model.nx_solver, model.nplist()), sdx_(model.nx_solver, model.nplist()),
      xB_(model.nJ * model.nx_solver), xQ_(model.nJ * model.nx_solver),
      xQB_(model.nplist()), xQBdot_(model.nplist()),
      dJydx_(model.nJ * model.nx_solver * model.nt(), 0.0) {
          /* maxSteps must be adapted if iterative linear solvers are used */
          if (solver.getLinearSolver() == LinearSolver::SPBCG) {
              max_steps_ = solver.getNewtonMaxSteps();
              numlinsteps_.resize(2 * max_steps_, 0);
          }
          /* Check for compatibility of options */
          if (solver.getSensitivityMethod() == SensitivityMethod::forward &&
              solver.getSensitivityMethodPreequilibration() == SensitivityMethod::adjoint &&
              solver.getSensitivityOrder() > SensitivityOrder::none)
              throw AmiException("Preequilibration using adjoint sensitivities "
                                 "is not compatible with using forward "
                                 "sensitivities during simulation");
      }

void SteadystateProblem::workSteadyStateProblem(Solver *solver, Model *model,
                                                int it) {

    /* process solver handling for pre- or postequilibration */
    if (it == -1) {
        /* solver was not run before, set up everything */
        model->initialize(x_, dx_, sx_, sdx_,
                          solver->getSensitivityOrder() >=
                              SensitivityOrder::first);
        t_ = model->t0();
        solver->setup(t_, model, x_, dx_, sx_, sdx_);
    } else {
        /* solver was run before, extract current state from solver */
        solver->writeSolution(&t_, x_, dx_, sx_, xQ_);
    }

    /* create a Newton solver object */
    auto newtonSolver = NewtonSolver::getSolver(&t_, &x_, *solver, model);

    /* Compute steady state and get the computation time */
    clock_t starttime = clock();
    findSteadyState(solver, newtonSolver.get(), model, it);
    cpu_time_ = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Check whether state sensis still need to be computed */
    if (getSensitivityFlag(model, solver, it, SteadyStateContext::newtonSensi))
    {
        try {
            /* this might still fail, if the Jacobian is singular and
             simulation did not find a steady state */
            newtonSolver->computeNewtonSensis(sx_);
        } catch (NewtonFailure const &) {
            /* No steady state could be inferred. Store simulation state */
            storeSimulationState(model, solver->getSensitivityOrder() >=
                                 SensitivityOrder::first);
            throw AmiException("Steady state sensitivity computation failed due "
                               "to unsuccessful factorization of RHS Jacobian");
        }
    }

    /* Get output of steady state solver, write it to x0 and reset time
     if necessary */
    storeSimulationState(model, getSensitivityFlag(model, solver, it,
                         SteadyStateContext::sensiStorage));
}

void SteadystateProblem::workSteadyStateBackwardProblem(Solver *solver,
                                                        Model *model,
                                                        const BackwardProblem *bwd) {
    /* initialize and check if there is something to be done */
    if (!initializeBackwardProblem(solver, model, bwd))
        return;

    /* Get the Newton solver */
    auto newtonSolver = NewtonSolver::getSolver(&t_, &x_, *solver, model);

    /* get the run time */
    clock_t starttime = clock();
    computeSteadyStateQuadrature(newtonSolver.get(), solver, model);
    cpu_timeB_ = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Finalize by setting adjoint state to zero (its steady state) */
    xB_.zero();
}

void SteadystateProblem::findSteadyState(Solver *solver,
                                         NewtonSolver *newtonSolver,
                                         Model *model, int it) {
    /* First, try to run the Newton solver */
    steady_state_status_.resize(3, SteadyStateStatus::not_run);
    findSteadyStateByNewtonsMethod(newtonSolver, model, false);

    /* Newton solver didn't work, so try to simulate to steady state */
    if (!checkSteadyStateSuccess())
        findSteadyStateBySimulation(solver, model, it);

    /* Simulation didn't work, retry the Newton solver from last sim state. */
    if (!checkSteadyStateSuccess())
        findSteadyStateByNewtonsMethod(newtonSolver, model, true);

    /* Nothing worked, throw an as informative error as possible */
    if (!checkSteadyStateSuccess())
        handleSteadyStateFailure(solver, model);
}

void SteadystateProblem::findSteadyStateByNewtonsMethod(NewtonSolver *newtonSolver,
                                                        Model *model,
                                                        bool newton_retry) {
    int ind = newton_retry ? 2 : 0;
    try {
        applyNewtonsMethod(model, newtonSolver, newton_retry);
        steady_state_status_[ind] = SteadyStateStatus::success;
    } catch (NewtonFailure const &ex) {
        /* nothing to be done */
        switch (ex.error_code) {
            case AMICI_TOO_MUCH_WORK:
                steady_state_status_[ind] =
                    SteadyStateStatus::failed_convergence;
                break;
            case AMICI_NO_STEADY_STATE:
                steady_state_status_[ind] =
                    SteadyStateStatus::failed_too_long_simulation;
                break;
            case AMICI_SINGULAR_JACOBIAN:
                steady_state_status_[ind] =
                    SteadyStateStatus::failed_factorization;
                break;
            case AMICI_DAMPING_FACTOR_ERROR:
                steady_state_status_[ind] = SteadyStateStatus::failed_damping;
                break;
            default:
                steady_state_status_[ind] = SteadyStateStatus::failed;
                break;
        }
    }

    /* copy number of linear steps used */
    if (max_steps_ > 0) {
        if (newton_retry) {
            std::copy_n(newtonSolver->getNumLinSteps().begin(),
                        max_steps_, &numlinsteps_.at(max_steps_));
        } else {
            std::copy_n(newtonSolver->getNumLinSteps().begin(),
                        max_steps_, numlinsteps_.begin());
        }
    }
}

void SteadystateProblem::findSteadyStateBySimulation(const Solver *solver,
                                                     Model *model,
                                                     int it) {
    /* set starting timepoint for the simulation solver */
    if (it < 1) /* No previous time point computed, set t = t0 */
        t_ = model->t0();
    else /* Carry on simulating from last point */
        t_ = model->getTimepoint(it - 1);

    try {
        if (it < 0) {
            /* Preequilibration? -> Create a new CVode object for sim */
            bool integrateSensis = getSensitivityFlag(model, solver, it,
                                   SteadyStateContext::solverCreation);
            auto newtonSimSolver = createSteadystateSimSolver(solver, model,
                                                              integrateSensis,
                                                              false);
            runSteadystateSimulation(newtonSimSolver.get(), model, false);
        } else {
            /* Solver was already created, use this one */
            runSteadystateSimulation(solver, model, false);
        }
        steady_state_status_[1] = SteadyStateStatus::success;
    } catch (NewtonFailure const &ex) {
        switch (ex.error_code) {
            case AMICI_TOO_MUCH_WORK:
                steady_state_status_[1] = SteadyStateStatus::failed_convergence;
                break;
            case AMICI_NO_STEADY_STATE:
                steady_state_status_[1] = SteadyStateStatus::failed_too_long_simulation;
                break;
            default:
                steady_state_status_[1] = SteadyStateStatus::failed;
        }
    } catch (AmiException const &) {
        steady_state_status_[1] = SteadyStateStatus::failed;
    }
}

bool SteadystateProblem::initializeBackwardProblem(Solver *solver,
                                                   Model *model,
                                                   const BackwardProblem *bwd) {
    if (bwd) {
        /* If preequilibration but not adjoint mode, there's nothing to do */
        if (solver->getSensitivityMethodPreequilibration() !=
            SensitivityMethod::adjoint)
            return false;

        /* If we need to reinitialize solver states, this won't work yet. */
        if (model->nx_reinit() > 0)
            throw NewtonFailure(AMICI_NOT_IMPLEMENTED,
                "Adjoint preequilibration with reinitialization of "
                "non-constant states is not yet implemented. Stopping.");

        /* If we have a backward problem, we're in preequilibration.
           Hence, quantities like t, x, and xB must be set. */
        solver->reInit(t_, x_, x_);
        solver->updateAndReinitStatesAndSensitivities(model);
        xB_.copy(bwd->getAdjointState());
    }

    /* Will need to write quadratures: set to 0 */
    xQ_.zero();
    xQB_.zero();
    xQBdot_.zero();

    return true;
}

void SteadystateProblem::computeSteadyStateQuadrature(NewtonSolver *newtonSolver,
                                                      const Solver *solver,
                                                      Model *model) {
    /* This routine computes the quadratures:
         xQB = Integral[ xB(x(t), t, p) * dxdot/dp(x(t), t, p) | dt ]
     As we're in steady state, we have x(t) = x_ss (x_steadystate), hence
         xQB = Integral[ xB(x_ss, t, p) | dt ] * dxdot/dp(x_ss, t, p)
     We therefore compute the integral over xB first and then do a
     matrix-vector multiplication */

    /* Try to compute the analytical solution for quadrature algebraically */
    getQuadratureByLinSolve(newtonSolver, model);

    /* Analytical solution didn't work, perform simulation instead */
    if (!hasQuadrature())
        getQuadratureBySimulation(solver, model);

    /* If analytic solution and integration did not work, throw an Exception */
    if (!hasQuadrature())
        throw AmiException("Steady state backward computation failed: Linear "
            "system could not be solved (possibly due to singular Jacobian), "
            "and numerical integration did not equilibrate within maxsteps");
}

void SteadystateProblem::getQuadratureByLinSolve(NewtonSolver *newtonSolver,
                                                 Model *model) {
    /* Computes the integral over the adjoint state xB:
     If the Jacobian has full rank, this has an analytical solution, since
     d/dt[ xB(t) ] = JB^T(x(t), p) xB(t) = JB^T(x_ss, p) xB(t)
     This linear ODE system with time-constant matrix has the solution
     xB(t) = exp( t * JB^T(x_ss, p) ) * xB(0)
     This integral xQ over xB is given as the solution of
     JB^T(x_ss, p) * xQ = xB(0)
     So we first try to solve the linear system, if possible. */

    /* copy content of xB into vector with integral */
    xQ_.copy(xB_);

    /* try to solve the linear system */
    try {
        /* compute integral over xB and write to xQ */
        newtonSolver->prepareLinearSystemB(0, -1);
        newtonSolver->solveLinearSystem(xQ_);
        /* Compute the quadrature as the inner product xQ * dxotdp */
        computeQBfromQ(model, xQ_, xQB_);
        /* set flag that quadratures is available (for processing in rdata) */
        hasQuadrature_ = true;
    } catch (NewtonFailure const &) {
        hasQuadrature_ = false;
    }
}

void SteadystateProblem::getQuadratureBySimulation(const Solver *solver,
                                                   Model *model) {
    /* If the Jacobian is singular, the integral over xB must be computed
       by usual integration over time, but  simplifications can be applied:
       x is not time dependent, no forward trajectory is needed. */

    /* set starting timepoint for the simulation solver */
    t_ = model->t0();
    /* xQ was written in getQuadratureByLinSolve() -> set to zero */
    xQ_.zero();

    /* create a new solver object */
    auto simSolver = createSteadystateSimSolver(solver, model, false, true);

    /* perform integration and quadrature */
    try {
        runSteadystateSimulation(simSolver.get(), model, true);
        hasQuadrature_ = true;
    } catch (NewtonFailure const &) {
        hasQuadrature_ = false;
    }
}

[[noreturn]] void SteadystateProblem::handleSteadyStateFailure(const Solver *solver,
                                                               Model *model) {
    /* No steady state could be inferred. Store simulation state */
    storeSimulationState(model, solver->getSensitivityOrder() >=
                         SensitivityOrder::first);

    /* Throw error message according to error codes */
    std::string errorString = "Steady state computation failed. "
                              "First run of Newton solver failed";
    writeErrorString(&errorString, steady_state_status_[0]);
    errorString.append(" Simulation to steady state failed");
    writeErrorString(&errorString, steady_state_status_[1]);
    errorString.append(" Second run of Newton solver failed");
    writeErrorString(&errorString, steady_state_status_[2]);

    throw AmiException(errorString.c_str());
}

void SteadystateProblem::writeErrorString(std::string *errorString,
                                          SteadyStateStatus status) const {
    /* write error message according to steady state status */
    switch (status) {
        case SteadyStateStatus::failed_too_long_simulation:
            (*errorString).append(": System could not be equilibrated via"
                                  " simulating to a late time point.");
            break;
        case SteadyStateStatus::failed_damping:
            (*errorString).append(": Damping factor reached lower bound.");
            break;
        case SteadyStateStatus::failed_factorization:
            (*errorString).append(": RHS could not be factorized.");
            break;
        case SteadyStateStatus::failed_convergence:
            (*errorString).append(": No convergence was achieved.");
            break;
        case SteadyStateStatus::failed:
            (*errorString).append(".");
            break;
        default:
            break;
    }
}

bool SteadystateProblem::getSensitivityFlag(const Model *model,
                                            const Solver *solver,
                                            int it, SteadyStateContext context) {
    /* We need to check whether we need to compute forward sensitivities.
       Depending on the situation (pre-/postequilibration) and the solver
       settings, the logic may be involved and is handled here.
       Most of these boolean operation could be simplified. However,
       clarity is more important than brevity. */

    /* Are we running in preequilibration (and hence create)? */
    bool preequilibration = (it == -1);

    /* Have we maybe already computed forward sensitivities? */
    bool forwardSensisAlreadyComputed =
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        steady_state_status_[1] == SteadyStateStatus::success &&
        model->getSteadyStateSensitivityMode() == SteadyStateSensitivityMode::simulationFSA;

    bool simulationStartedInSteadystate =
        steady_state_status_[0] == SteadyStateStatus::success &&
        numsteps_[0] == 0;

    /* Do we need forward sensis for postequilibration? */
    bool needForwardSensisPosteq = !preequilibration &&
        !forwardSensisAlreadyComputed &&
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        solver->getSensitivityMethod() == SensitivityMethod::forward;

    /* Do we need forward sensis for preequilibration? */
    bool needForwardSensisPreeq = preequilibration &&
        !forwardSensisAlreadyComputed &&
        solver->getSensitivityMethodPreequilibration() == SensitivityMethod::forward &&
        solver->getSensitivityOrder() >= SensitivityOrder::first;

    /* Do we need to do the linear system solve to get forward sensitivities? */
    bool needForwardSensisNewton =
        (needForwardSensisPreeq || needForwardSensisPosteq) &&
        !simulationStartedInSteadystate;

    /* When we're creating a new solver object */
    bool needForwardSensiAtCreation = needForwardSensisPreeq &&
        model->getSteadyStateSensitivityMode() == SteadyStateSensitivityMode::simulationFSA;

    /* Check if we need to store sensis */
    switch (context) {
        case SteadyStateContext::newtonSensi:
            return needForwardSensisNewton;

        case SteadyStateContext::sensiStorage:
            return needForwardSensisNewton ||
                forwardSensisAlreadyComputed ||
                simulationStartedInSteadystate;

        case SteadyStateContext::solverCreation:
            return needForwardSensiAtCreation;

        default:
            throw AmiException("Requested invalid context in sensitivity "
                               "processing during steady state computation");
    }
}

realtype SteadystateProblem::getWrmsNorm(const AmiVector &x,
                                         const AmiVector &xdot,
                                         realtype atol,
                                         realtype rtol,
                                         AmiVector &ewt) const {
    /* Depending on what convergence we want to check (xdot, sxdot, xQBdot)
       we need to pass ewt[QB], as xdot and xQBdot have different sizes */
    N_VAbs(const_cast<N_Vector>(x.getNVector()), ewt.getNVector());
    N_VScale(rtol, ewt.getNVector(), ewt.getNVector());
    N_VAddConst(ewt.getNVector(), atol, ewt.getNVector());
    N_VInv(ewt.getNVector(), ewt.getNVector());
    return N_VWrmsNorm(const_cast<N_Vector>(xdot.getNVector()),
                       ewt.getNVector());
}

bool SteadystateProblem::checkConvergence(const Solver *solver,
                                          Model *model,
                                          SensitivityMethod checkSensitivities) {
    if (checkSensitivities == SensitivityMethod::adjoint) {
        /* In the adjoint case, only xQB contributes to the gradient, the exact
           steadystate is less important, as xB = xQdot may even not converge
           to zero at all. So we need xQBdot, hence compute xQB first. */
        computeQBfromQ(model, xQ_, xQB_);
        computeQBfromQ(model, xB_, xQBdot_);
        wrms_ = getWrmsNorm(xQB_, xQBdot_, solver->getAbsoluteToleranceQuadratures(),
                           solver->getRelativeToleranceQuadratures(), ewtQB_);
    } else {
        /* If we're doing a forward simulation (with or without sensitivities:
           Get RHS and compute weighted error norm */
        model->fxdot(t_, x_, dx_, xdot_);
        wrms_ = getWrmsNorm(x_, xdot_, solver->getAbsoluteToleranceSteadyState(),
                           solver->getRelativeToleranceSteadyState(), ewt_);
    }
    bool converged = wrms_ < RCONST(1.0);

    if (checkSensitivities != SensitivityMethod::forward)
        return converged;

    /* Forward sensitivities: Compute weighted error norm for their RHS */
    for (int ip = 0; ip < model->nplist(); ++ip) {
        if (converged) {
            sx_ = solver->getStateSensitivity(t_);
            model->fsxdot(t_, x_, dx_, ip, sx_[ip], dx_, xdot_);
            wrms_ = getWrmsNorm(
                x_, xdot_, solver->getAbsoluteToleranceSteadyStateSensi(),
                solver->getRelativeToleranceSteadyStateSensi(), ewt_);
            converged = wrms_ < RCONST(1.0);
        }
    }
    return converged;
}

bool SteadystateProblem::checkSteadyStateSuccess() const {
    /* Did one of the attempts yield s steady state? */
    if (std::any_of(steady_state_status_.begin(), steady_state_status_.end(),
                    [](SteadyStateStatus status)
                    {return status == SteadyStateStatus::success;})) {
        return true;
    } else {
        return false;
    }
}

void SteadystateProblem::applyNewtonsMethod(Model *model,
                                            NewtonSolver *newtonSolver,
                                            bool newton_retry) {
    int i_newtonstep = 0;
    int ix = 0;
    double gamma = 1.0;
    bool compNewStep = true;

    if (model->nx_solver == 0)
        return;

    /* initialize output of linear solver for Newton step */
    delta_.zero();

    model->fxdot(t_, x_, dx_, xdot_);

    /* Check for relative error, but make sure not to divide by 0!
        Ensure positivity of the state */
    x_newton_ = x_;
    x_old_ = x_;
    xdot_old_ = xdot_;

    wrms_ = getWrmsNorm(x_newton_, xdot_, newtonSolver->atol_,
                       newtonSolver->rtol_, ewt_);
    bool converged = wrms_ < RCONST(1.0);
    while (!converged && i_newtonstep < newtonSolver->max_steps) {

        /* If Newton steps are necessary, compute the initial search direction */
        if (compNewStep) {
            try {
                delta_ = xdot_;
                newtonSolver->getStep(newton_retry ? 2 : 1, i_newtonstep, delta_);
            } catch (NewtonFailure const &) {
                numsteps_.at(newton_retry ? 2 : 0) = i_newtonstep;
                throw;
            }
        }

        /* Try a full, undamped Newton step */
        N_VLinearSum(1.0, x_old_.getNVector(), gamma, delta_.getNVector(),
                     x_.getNVector());

        /* Compute new xdot and residuals */
        model->fxdot(t_, x_, dx_, xdot_);
        realtype wrms_tmp = getWrmsNorm(x_newton_, xdot_, newtonSolver->atol_,
                                        newtonSolver->rtol_, ewt_);

        if (wrms_tmp < wrms_) {
            /* If new residuals are smaller than old ones, update state */
            wrms_ = wrms_tmp;
            x_old_ = x_;
            xdot_old_ = xdot_;
            /* New linear solve due to new state */
            compNewStep = true;
            /* Check residuals vs tolerances */
            converged = wrms_ < RCONST(1.0);

            if (converged) {
                /* Ensure positivity of the found state and recheck if
                   the convergence still holds */
                bool recheck_convergence = false;
                for (ix = 0; ix < model->nx_solver; ix++) {
                    if (x_[ix] < 0.0) {
                        x_[ix] = 0.0;
                        recheck_convergence = true;
                    }
                }
                if (recheck_convergence) {
                  model->fxdot(t_, x_, dx_, xdot_);
                  wrms_ = getWrmsNorm(x_newton_, xdot_, newtonSolver->atol_,
                                     newtonSolver->rtol_, ewt_);
                  converged = wrms_ < RCONST(1.0);
                }
            } else if (newtonSolver->damping_factor_mode_==NewtonDampingFactorMode::on) {
                /* increase dampening factor (superfluous, if converged) */
                gamma = fmin(1.0, 2.0 * gamma);
            }
        } else if (newtonSolver->damping_factor_mode_==NewtonDampingFactorMode::on) {
            /* Reduce dampening factor and raise an error when becomes too small */
            gamma = gamma / 4.0;
            if (gamma < newtonSolver->damping_factor_lower_bound)
              throw NewtonFailure(AMICI_DAMPING_FACTOR_ERROR,
                                  "Newton solver failed: the damping factor "
                                  "reached its lower bound");

            /* No new linear solve, only try new dampening */
            compNewStep = false;
        }
        /* increase step counter */
        i_newtonstep++;
    }

    /* Set return values */
    numsteps_.at(newton_retry ? 2 : 0) = i_newtonstep;
    if (!converged)
        throw NewtonFailure(AMICI_TOO_MUCH_WORK, "applyNewtonsMethod");
}

void SteadystateProblem::runSteadystateSimulation(const Solver *solver,
                                                  Model *model,
                                                  bool backward)
{
    if (model->nx_solver == 0)
        return;
    /* Loop over steps and check for convergence
       NB: This function is used for forward and backward simulation, and may
       be called by workSteadyStateProblem and workSteadyStateBackwardProblem.
       Whether we simulate forward or backward in time is reflected by the
       flag "backward". */

    /* Do we also have to check for convergence of sensitivities? */
    SensitivityMethod sensitivityFlag = SensitivityMethod::none;
    if (solver->getSensitivityOrder() > SensitivityOrder::none &&
        solver->getSensitivityMethod() > SensitivityMethod::none)
        sensitivityFlag = SensitivityMethod::forward;
    /* If flag for forward sensitivity computation by simulation is not set,
     disable forward sensitivity integration. Sensitivities will be computed
     by newtonSolver->computeNewtonSensis then */
    if (model->getSteadyStateSensitivityMode() == SteadyStateSensitivityMode::newtonOnly) {
        solver->switchForwardSensisOff();
        sensitivityFlag = SensitivityMethod::none;
    }
    if (backward)
        sensitivityFlag = SensitivityMethod::adjoint;

    bool converged = checkConvergence(solver, model, sensitivityFlag);
    int sim_steps = 0;

    while (!converged) {
        /* One step of ODE integration
         reason for tout specification:
         max with 1 ensures correct direction (any positive value would do)
         multiplication with 10 ensures nonzero difference and should ensure
         stable computation value is not important for AMICI_ONE_STEP mode,
         only direction w.r.t. current t
         */
        solver->step(std::max(t_, 1.0) * 10);
        if (backward) {
            solver->writeSolution(&t_, xB_, dx_, sx_, xQ_);
        } else {
            solver->writeSolution(&t_, x_, dx_, sx_, xQ_);
        }

        /* Check for convergence */
        converged = checkConvergence(solver, model, sensitivityFlag);
        /* increase counter, check for maxsteps */
        sim_steps++;
        if (sim_steps >= solver->getMaxSteps() && !converged) {
            numsteps_.at(1) = sim_steps;
            throw NewtonFailure(AMICI_TOO_MUCH_WORK,
                                "exceeded maximum number of steps");
        }
        if (t_ >= 1e200 && !converged) {
            numsteps_.at(1) = sim_steps;
            throw NewtonFailure(AMICI_NO_STEADY_STATE, "simulated to late time"
                                " point without convergence of RHS");
        }
    }

    /* store information about steps and sensitivities, if necessary */
    if (backward) {
        numstepsB_ = sim_steps;
    } else {
        numsteps_.at(1) = sim_steps;
        if (solver->getSensitivityOrder() > SensitivityOrder::none &&
            model->getSteadyStateSensitivityMode() ==
            SteadyStateSensitivityMode::simulationFSA)
            sx_ = solver->getStateSensitivity(t_);
    }
}

std::unique_ptr<Solver> SteadystateProblem::createSteadystateSimSolver(
        const Solver *solver, Model *model, bool forwardSensis, bool backward) const
{
    /* Create new CVode solver object */
    auto sim_solver = std::unique_ptr<Solver>(solver->clone());

    switch (solver->getLinearSolver()) {
        case LinearSolver::dense:
            break;
        case LinearSolver::KLU:
            break;
        default:
            throw NewtonFailure(AMICI_NOT_IMPLEMENTED,
                                "invalid solver for steadystate simulation");
    }
    /* do we need sensitivities? */
    if (forwardSensis) {
        /* need forward to compute sx0 */
        sim_solver->setSensitivityMethod(SensitivityMethod::forward);
    } else {
        sim_solver->setSensitivityMethod(SensitivityMethod::none);
        sim_solver->setSensitivityOrder(SensitivityOrder::none);
    }
    /* use x and sx as dummies for dx and sdx
     (they wont get touched in a CVodeSolver) */
    sim_solver->setup(model->t0(), model, x_, x_, sx_, sx_);
    if (backward) {
        sim_solver->setup(model->t0(), model, xB_, xB_, sx_, sx_);
        sim_solver->setupSteadystate(model->t0(), model, x_, x_, xB_, xB_, xQ_);
    } else {
        sim_solver->setup(model->t0(), model, x_, x_, sx_, sx_);
    }

    return sim_solver;
}

void SteadystateProblem::computeQBfromQ(Model *model, const AmiVector &yQ,
                                        AmiVector &yQB) const {
    /* Compute the quadrature as the inner product: yQB = dxotdp * yQ */

    /* set to zero first, as multiplication adds to existing value */
    yQB.zero();
    /* multiply */
    if (model->pythonGenerated) {
        /* fill dxdotdp with current values */
        const auto& plist = model->getParameterList();
        model->fdxdotdp(t_, x_, x_);
        model->get_dxdotdp_full().multiply(yQB.getNVector(), yQ.getNVector(),
                                           plist, true);
    } else {
        for (int ip=0; ip<model->nplist(); ++ip)
            yQB[ip] = N_VDotProd(
                const_cast<N_Vector>(yQ.getNVector()),
                const_cast<N_Vector>(model->get_dxdotdp().getNVector(ip)));
    }
}

void SteadystateProblem::getAdjointUpdates(Model &model,
                                           const ExpData &edata) {
    xB_.zero();
    for (int it=0; it < model.nt(); it++) {
        if (std::isinf(model.getTimepoint(it))) {
            model.getAdjointStateObservableUpdate(
                slice(dJydx_, it, model.nx_solver * model.nJ), it, x_, edata);
            for (int ix = 0; ix < model.nxtrue_solver; ix++)
                xB_[ix] += dJydx_[ix + it * model.nx_solver];
        }
    }
}

void SteadystateProblem::storeSimulationState(Model *model, bool storesensi) {
    state_.t = INFINITY;
    state_.x = x_;
    state_.dx = xdot_;
    if (storesensi)
        state_.sx = sx_;
    state_.state = model->getModelState();
}

} // namespace amici
