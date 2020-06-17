#include "amici/steadystateproblem.h"
#include "amici/defines.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/solver_cvodes.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
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
    : delta(model.nx_solver), ewt(model.nx_solver),
      rel_x_newton(model.nx_solver), x_newton(model.nx_solver),
      x(model.nx_solver), x_old(model.nx_solver), dx(model.nx_solver),
      xdot(model.nx_solver), xdot_old(model.nx_solver),
      sx(model.nx_solver, model.nplist()), sdx(model.nx_solver, model.nplist()),
      dJydx(model.nJ * model.nx_solver * model.nt(), 0.0), numsteps(3, 0) {
          /* maxSteps must be adapted if iterative linear solvers are used */
          if (solver.getLinearSolver() == LinearSolver::SPBCG) {
              maxSteps = solver.getNewtonMaxSteps();
              numlinsteps.resize(2 * maxSteps, 0);
          }
      }

void SteadystateProblem::workSteadyStateProblem(Solver *solver, Model *model,
                                                int it) {

    /* process solver handling for pre- or postequilibration */
    if (it == -1) {
        /* solver was not run before, set up everything */
        model->initialize(x, dx, sx, sdx,
                          solver->getSensitivityOrder() >=
                              SensitivityOrder::first);
        t = model->t0();
        solver->setup(t, model, x, dx, sx, sdx);
    } else {
        /* Are we computing adjoint sensitivities? That's not yet supoorted
           with steady state sensitivity analysis */
        if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
            solver->getSensitivityMethod() == SensitivityMethod::adjoint)
            throw AmiException("Steady state sensitivity computation together "
                               "with adjoint sensitivity analysis is currently "
                               "not supported.");
        /* solver was run before, extract current state from solver */
        solver->writeSolution(&t, x, dx, sx);
    }

    /* create a Newton solver obejct */
    auto newtonSolver = NewtonSolver::getSolver(&t, &x, *solver, model);

    /* Compute steady state and get the computation time */
    clock_t starttime = clock();
    findSteadyState(solver, newtonSolver.get(), model, it);
    cpu_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Check whether state sensis still need to be computed */
    if (processSensitivityLogic(model)) {
        try {
            /* this might still fail, if the Jacobian is singular and
             simulation did not find a steady state */
            newtonSolver->computeNewtonSensis(sx);
        } catch (NewtonFailure const &) {
            /* No steady state could be inferred. Store simulation state */
            storeSimulationState(model, solver->getSensitivityOrder() >=
                                 SensitivityOrder::first);
            throw AmiException("Steady state sensitvitiy computation failed due "
                               "to unsuccessful factorization of RHS Jacobian");
        }
    }

    /* Get output of steady state solver, write it to x0 and reset time
     if necessary */
    storeSimulationState(model, solver->getSensitivityOrder() >=
                                    SensitivityOrder::first);
}

void SteadystateProblem::findSteadyState(Solver *solver,
                                         NewtonSolver *newtonSolver,
                                         Model *model, int it) {
    /* First, try to run the Newton solver */
    findSteadyStateByNewtonsMethod(newtonSolver, model, false);

    /* Newton solver didn't work, so try to simulate to steady state */
    if (!checkSteadyStateStatus())
        findSteadyStateBySimulation(solver, model, it);

    /* Simulation didn't work, retry the Newton solver from last sim state. */
    if (!checkSteadyStateStatus())
        findSteadyStateByNewtonsMethod(newtonSolver, model, true);

    /* Nothing worked, throw an as informative error as possible */
    if (!checkSteadyStateStatus())
        handleSteadyStateComputationFailure(model, solver);
}

void SteadystateProblem::findSteadyStateByNewtonsMethod(NewtonSolver *newtonSolver,
                                                        Model *model,
                                                        bool newton_retry) {
    try {
        applyNewtonsMethod(model, newtonSolver, newton_retry);
        if (newton_retry)
            steady_state_status[2] = SteadyStateStatus::success;
        else
            steady_state_status[0] = SteadyStateStatus::success;
    } catch (NewtonFailure const &ex) {
        /* nothing to be done */
        switch (ex.error_code) {
            case AMICI_TOO_MUCH_WORK:
                steady_state_status[0] = SteadyStateStatus::failed_convergence;
                break;
            case AMICI_SINGULAR_JACOBIAN:
                steady_state_status[0] = SteadyStateStatus::failed_factorization;
                break;
            case AMICI_DAMPING_FACTOR_ERROR:
                steady_state_status[0] = SteadyStateStatus::failed_damping;
                break;
            default:
                steady_state_status[0] = SteadyStateStatus::failed;
    }

    /* copy number of linear steps used */
    if (maxSteps > 0)
        std::copy_n(newtonSolver->getNumLinSteps().begin(),
                    maxSteps, numlinsteps.begin());
}

void SteadystateProblem::findSteadyStateByNewtonSimulation(Solver *solver,
                                                           Model *model,
                                                           int it) {
    /* set starting timepoint for the simulation solver */
    if (it < 1) /* No previous time point computed, set t = t0 */
        t = model->t0();
    else /* Carry on simulating from last point */
        t = model->getTimepoint(it - 1);

    try {
        if (it < 0) {
            /* Preequilibration? -> Create a new CVode object for sim */
            auto newtonSimSolver = createSteadystateSimSolver(solver, model);
            getSteadystateSimulation(newtonSimSolver.get(), model);
        } else {
            /* Solver was already created, use this one */
            getSteadystateSimulation(solver, model);
        }
        steady_state_status[1] = SteadyStateStatus::success;
    } catch (NewtonFailure const &ex) {
        if (ex.error_code == AMICI_TOO_MUCH_WORK) {
            steady_state_status[1] = SteadyStateStatus::failed_convergence;
        } else {
            steady_state_status[1] = SteadyStateStatus::failed;
        }
    } catch (AmiException const &) {
        steady_state_status[1] = SteadyStateStatus::failed;
    }
}

void SteadystateProblem::handleSteadyStateComputationFailure(Solver *solver,
                                                             Model *model) {
    /* No steady state could be inferred. Store simulation state */
    storeSimulationState(model, solver->getSensitivityOrder() >=
                         SensitivityOrder::first);

    /* Throw error message according to error codes */
    string error_string = "Steady state computation failed. "
                          "First run of Newton solver failed";
    error_string = write_error_string(error_string, steady_state_status[0]);
    error_string.apppend(" Simulation to steady state failed");
    error_string = write_error_string(error_string, steady_state_status[1]);
    error_string.apppend(" Second run of Newton solver failed");
    error_string = write_error_string(error_string, steady_state_status[2]);

    throw AmiException(error_string);
}

std::string SteadystateProblem::write_error_string(std::string error_string,
                                                   SteadyStateStatus status) {
    /* write error message according to steady state status */
    switch (status) {
        case SteadyStateStatus::failed_damping:
            error_string.append(": Damping factor reached lower bound.");
        case SteadyStateStatus::failed_factorization:
            error_string.append(": RHS could not be factorized.");
        case SteadyStateStatus::failed_factorization:
            error_string.append(": No convergence was achieved.");
        case SteadyStateStatus::failed:
            error_string.append(".");
    }
    return error_string;
}


bool SteadystateProblem::processSensitivityLogic(Model *model) {
    /* NB: Currently, this logic processing is still "simple".
           However, it will become more involved once adjoint
           sensitivities are implemented */

    /* We want to solve the linear system if newtonOnly was used... */
    bool needStateSenis = model->getSteadyStateSensitivityMode() ==
        SteadyStateSensitivityMode::newtonOnly;
    /* ... or if Newton's method found the steady state */
    needStateSenis = needStateSenis ||
        steady_state_status[0] == SteadyStateStatus::success ||
        steady_state_status[2] == SteadyStateStatus::success;
    return needStateSenis;
}

realtype SteadystateProblem::getWrmsNorm(const AmiVector &x,
                                         const AmiVector &xdot,
                                         realtype atol,
                                         realtype rtol
                                         ) {
    N_VAbs(x.getNVector(), ewt.getNVector());
    N_VScale(rtol, ewt.getNVector(), ewt.getNVector());
    N_VAddConst(ewt.getNVector(), atol, ewt.getNVector());
    N_VInv(ewt.getNVector(), ewt.getNVector());
    return N_VWrmsNorm(xdot.getNVector(), ewt.getNVector());
}

bool SteadystateProblem::checkConvergence(const Solver *solver, Model *model) {
    /* get RHS and compute weighted error norm */
    model->fxdot(t, x, dx, xdot);
    wrms = getWrmsNorm(x, xdot, solver->getAbsoluteToleranceSteadyState(),
                       solver->getRelativeToleranceSteadyState());
    bool converged = wrms < RCONST(1.0);

    /* If we also integrate forward sensis, we need to also check those:
       Check if: sensis enabled && steadyStateSensiMode == simulation
     */
    bool checkForwardSensis =
        solver->getSensitivityOrder() > SensitivityOrder::none &&
        model->getSteadyStateSensitivityMode() ==
            SteadyStateSensitivityMode::simulationFSA;

    /* get RHS of sensitivities and compute weighted error norm */
    if (checkForwardSensis) {
        for (int ip = 0; ip < model->nplist(); ++ip) {
            if (converged) {
                sx = solver->getStateSensitivity(t);
                model->fsxdot(t, x, dx, ip, sx[ip], dx, xdot);
                wrms = getWrmsNorm(
                    x, xdot, solver->getAbsoluteToleranceSteadyStateSensi(),
                    solver->getRelativeToleranceSteadyStateSensi());
                converged = wrms < RCONST(1.0);
            }
        }
    }
    return converged;
}

bool checkSteadyStateSuccess() {
    /* Did one of the attempts yield s steady state? */
    if (steady_state_status[0] == SteadyStateStatus::success ||
        steady_state_status[1] == SteadyStateStatus::success ||
        steady_state_status[2] == SteadyStateStatus::success) {
        return true;
    } else {
        return false;
    }
};

void SteadystateProblem::applyNewtonsMethod(Model *model,
                                            NewtonSolver *newtonSolver,
                                            bool newton_retry) {
    int i_newtonstep = 0;
    int ix = 0;
    double gamma = 1.0;
    bool compNewStep = true;

    /* initialize output of linear solver for Newton step */
    delta.reset();

    model->fxdot(t, x, dx,xdot);

    /* Check for relative error, but make sure not to divide by 0!
        Ensure positivity of the state */
    x_newton = x;
    x_old = x;
    xdot_old = xdot;

    wrms = getWrmsNorm(x_newton, xdot, newtonSolver->atol, newtonSolver->rtol);
    bool converged = wrms < RCONST(1.0);
    while (!converged && i_newtonstep < newtonSolver->maxsteps) {

        /* If Newton steps are necessary, compute the inital search direction */
        if (compNewStep) {
            try {
                delta = xdot;
                newtonSolver->getStep(newton_retry ? 2 : 1, i_newtonstep, delta);
            } catch (NewtonFailure const &) {
                numsteps.at(newton_retry ? 2 : 0) = i_newtonstep;
                throw;
            }
        }

        /* Try a full, undamped Newton step */
        N_VLinearSum(1.0, x_old.getNVector(), gamma, delta.getNVector(),
                     x.getNVector());

        /* Compute new xdot and residuals */
        model->fxdot(t, x, dx, xdot);
        realtype wrms_tmp = getWrmsNorm(x_newton, xdot, newtonSolver->atol,
                                        newtonSolver->rtol);

        if (wrms_tmp < wrms) {
            /* If new residuals are smaller than old ones, update state */
            wrms = wrms_tmp;
            x_old = x;
            xdot_old = xdot;
            /* New linear solve due to new state */
            compNewStep = true;
            /* Check residuals vs tolerances */
            converged = wrms < RCONST(1.0);

            if (converged) {
                /* Ensure positivity of the found state and recheck if
                   the convergence still holds */
                bool recheck_convergence = false;
                for (ix = 0; ix < model->nx_solver; ix++) {
                    if (x[ix] < 0.0) {
                        x[ix] = 0.0;
                        recheck_convergence = true;
                    }
                }
                if (recheck_convergence) {
                  model->fxdot(t, x, dx, xdot);
                  wrms = getWrmsNorm(x_newton, xdot, newtonSolver->atol, newtonSolver->rtol);
                  converged = wrms < RCONST(1.0);
                }
            } else if (newtonSolver->dampingFactorMode==NewtonDampingFactorMode::on) {
                /* increase dampening factor (superfluous, if converged) */
                gamma = fmin(1.0, 2.0 * gamma);
            }
        } else if (newtonSolver->dampingFactorMode==NewtonDampingFactorMode::on) {
            /* Reduce dampening factor and raise an error when becomes too small */
            gamma = gamma / 4.0;
            if (gamma < newtonSolver->dampingFactorLowerBound)
              throw NewtonFailure(AMICI_CONV_FAILURE,
                                  "Newton solver failed: the damping factor "
                                  "reached its lower bound");

            /* No new linear solve, only try new dampening */
            compNewStep = false;
        }
        /* increase step counter */
        i_newtonstep++;
    }

    /* Set return values */
    numsteps.at(newton_retry ? 2 : 0) = i_newtonstep;
    if (!converged)
        throw NewtonFailure(AMICI_TOO_MUCH_WORK, "applyNewtonsMethod");
}

void SteadystateProblem::getSteadystateSimulation(Solver *solver,
                                                  Model *model)
{
    /* Loop over steps and check for convergence */
    bool converged = checkConvergence(solver, model);
    int sim_steps = 0;
    /* If flag for forward sensitivity computation by simulation is not set,
     disable forward sensitivity integration. Sensitivities will be combputed
     by newonSolver->computeNewtonSensis then */
    if (model->getSteadyStateSensitivityMode() ==
        SteadyStateSensitivityMode::newtonOnly)
        solver->switchForwardSensisOff();

    while (!converged) {
        /* One step of ODE integration
         reason for tout specification:
         max with 1 ensures correct direction (any positive value would do)
         multiplication with 10 ensures nonzero difference and should ensure
         stable computation value is not important for AMICI_ONE_STEP mode,
         only direction w.r.t. current t
         */
        solver->step(std::max(t, 1.0) * 10);
        solver->writeSolution(&t, x, dx, sx);

        /* Check for convergence */
        converged = checkConvergence(solver, model);
        /* increase counter, check for maxsteps */
        sim_steps++;
        if (sim_steps >= solver->getMaxSteps() && !converged) {
            numsteps.at(1) = sim_steps;
            throw NewtonFailure(AMICI_TOO_MUCH_WORK,
                                "exceeded maximum number of steps");
        }
        if (t >= 1e100 && !converged) {
            numsteps.at(1) = sim_steps;
            throw NewtonFailure(AMICI_TOO_MUCH_WORK,
                                "simulated beyond t=1e100 without convergence");
        }
    }
    numsteps.at(1) = sim_steps;
    if (solver->getSensitivityOrder() > SensitivityOrder::none &&
        model->getSteadyStateSensitivityMode() ==
        SteadyStateSensitivityMode::simulationFSA)
        sx = solver->getStateSensitivity(t);
}

std::unique_ptr<Solver> SteadystateProblem::createSteadystateSimSolver(
        const Solver *solver, Model *model) const
{
    /* Create new CVode solver object */

    auto newton_solver = std::unique_ptr<Solver>(solver->clone());

    switch (solver->getLinearSolver()) {
    case LinearSolver::dense:
    case LinearSolver::KLU:
    case LinearSolver::SuperLUMT:
        break;
    default:
        throw NewtonFailure(AMICI_NOT_IMPLEMENTED,
                            "invalid solver for steadystate simulation");
    }
    if (solver->getSensitivityMethod() != SensitivityMethod::none &&
        model->getSteadyStateSensitivityMode() ==
            SteadyStateSensitivityMode::simulationFSA)
        newton_solver->setSensitivityMethod(SensitivityMethod::forward);
    // need forward to compute sx0
    else
        newton_solver->setSensitivityMethod(SensitivityMethod::none);

    /* use x and sx as dummies for dx and sdx
     (they wont get touched in a CVodeSolver) */
    newton_solver->setup(model->t0(), model, x, x, sx, sx);

    return newton_solver;
}

void SteadystateProblem::getAdjointUpdates(Model &model,
                                           const ExpData &edata) {
    for (int it=0; it < model.nt(); it++) {
        if (std::isinf(model.getTimepoint(it))) {
            model.getAdjointStateObservableUpdate(
                slice(dJydx, it, model.nx_solver * model.nJ), it, x, edata);
        }
    }
}

void SteadystateProblem::storeSimulationState(Model *model, bool storesensi) {
    state.t = INFINITY;
    state.x = x;
    state.dx = xdot;
    if (storesensi)
        state.sx = sx;
    state.state = model->getModelState();
}

} // namespace amici
