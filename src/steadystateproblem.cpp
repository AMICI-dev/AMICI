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
    : delta(model.nx_solver), ewt(model.nx_solver),
      rel_x_newton(model.nx_solver), x_newton(model.nx_solver),
      x(model.nx_solver), x_old(model.nx_solver), dx(model.nx_solver),
      xdot(model.nx_solver), xdot_old(model.nx_solver),
      sx(model.nx_solver, model.nplist()), sdx(model.nx_solver, model.nplist()),
      xB(model.nJ * model.nx_solver), xQB(model.nplist()),
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
    if (getSensitivityFlag(model, solver, it, SteadyStateContext::newton))
    {
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
    storeSimulationState(model, getSensitivityFlag(model, solver, it,
                                                   SteadyStateContext::storage));
}

void SteadystateProblem::workSteadyStateBackwardProblem(Solver *solver,
                                                        Model *model,
                                                        BackwardProblem *bwd) {
    /* initialize and check if there is something to be done */
    bool computeBackward = initializeBackwardProblem(solver, model, bwd);
    if (!computeBackward)
        return;

    /* Get the Newton solver */
    auto newtonSolver = NewtonSolver::getSolver(&t, &x, *solver, model);

    /* get the run time */
    clock_t starttime = clock();
    computeSteadyStateQuadrature(newtonSolver.get(), model);
    cpu_timeB = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Finalize by setting addjoint state to zero (its steady state) */
    xB.reset();
}

void SteadystateProblem::findSteadyState(Solver *solver,
                                         NewtonSolver *newtonSolver,
                                         Model *model, int it) {
    /* First, try to run the Newton solver */
    steady_state_status.resize(3, SteadyStateStatus::not_run);
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
        steady_state_status[ind] = SteadyStateStatus::success;
    } catch (NewtonFailure const &ex) {
        /* nothing to be done */
        switch (ex.error_code) {
            case AMICI_TOO_MUCH_WORK:
                steady_state_status[ind] =
                    SteadyStateStatus::failed_convergence;
                break;
            case AMICI_NO_STEADY_STATE:
                steady_state_status[ind] =
                    SteadyStateStatus::failed_too_long_simulation;
                break;
            case AMICI_SINGULAR_JACOBIAN:
                steady_state_status[ind] =
                    SteadyStateStatus::failed_factorization;
                break;
            case AMICI_DAMPING_FACTOR_ERROR:
                steady_state_status[ind] = SteadyStateStatus::failed_damping;
                break;
            default:
                steady_state_status[ind] = SteadyStateStatus::failed;
                break;
        }
    }

    /* copy number of linear steps used */
    if (maxSteps > 0) {
        if (newton_retry) {
            std::copy_n(newtonSolver->getNumLinSteps().begin(),
                        maxSteps, &numlinsteps.at(maxSteps));
        } else {
            std::copy_n(newtonSolver->getNumLinSteps().begin(),
                        maxSteps, numlinsteps.begin());
        }
    }
}

void SteadystateProblem::findSteadyStateBySimulation(Solver *solver,
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
            bool integrateSensis =
                getSensitivityFlag(model, solver, it,
                                   SteadyStateContext::integration);
            auto newtonSimSolver = createSteadystateSimSolver(solver, model,
                                                              integrateSensis);
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

bool SteadystateProblem::initializeBackwardProblem(Solver *solver,
                                                   Model *model,
                                                   BackwardProblem *bwd) {
    if (bwd) {
        /* If preequilibration but not adjoint mode, there's nothing to do */
        if (solver->getSensitivityMethodPreequilibration() !=
            SensitivityMethod::adjoint) {
            return false;
        }

        /* If we have a backward problem, we're in preequilibration.
           Hence, quantities like t, x, and xB must be set. */
        solver->updateAndReinitStatesAndSensitivities(model);
        xB.copy(bwd->getAdjointState());
    } else {
        xB.reset();
    }

    /* Will need to write quadratures: resize */
    xQB.reset();
    return true;
}

void SteadystateProblem::computeSteadyStateQuadrature(NewtonSolver *newtonSolver,
                                                      Model *model) {
    /* compute the integral over the adjoint state xBintegral */
    try {
        newtonSolver->prepareLinearSystemB(0, -1);
        newtonSolver->solveLinearSystem(xB);
    } catch (NewtonFailure const &ex) {
        if (ex.error_code == AMICI_SINGULAR_JACOBIAN)
            throw NewtonFailure(ex.error_code, "Steady state backward "
                                "computation failed due to unsuccessful "
                                "factorization of RHS Jacobian. ");
        throw NewtonFailure(ex.error_code, "Steady state backward "
                            "computation failed. ");
    }

    /* Compute the quadrature as the inner product xBintegral * dxotdp */
    if (model->pythonGenerated) {
        const auto& plist = model->getParameterList();
        if (model->ndxdotdp_explicit > 0)
            model->dxdotdp_explicit.multiply(xQB.getNVector(),
                                             xB.getNVector(), plist, true);
        if (model->ndxdotdp_implicit > 0)
            model->dxdotdp_implicit.multiply(xQB.getNVector(),
                                             xB.getNVector(), plist, true);
    } else {
        for (int ip=0; ip<model->nplist(); ++ip)
            xQB[ip] = N_VDotProd(xB.getNVector(),
                                 model->dxdotdp.getNVector(ip));
    }
    hasQuadrature = true;
}

[[noreturn]] void SteadystateProblem::handleSteadyStateFailure(const Solver *solver,
                                                               Model *model) {
    /* No steady state could be inferred. Store simulation state */
    storeSimulationState(model, solver->getSensitivityOrder() >=
                         SensitivityOrder::first);

    /* Throw error message according to error codes */
    std::string errorString = "Steady state computation failed. "
                              "First run of Newton solver failed";
    writeErrorString(&errorString, steady_state_status[0]);
    errorString.append(" Simulation to steady state failed");
    writeErrorString(&errorString, steady_state_status[1]);
    errorString.append(" Second run of Newton solver failed");
    writeErrorString(&errorString, steady_state_status[2]);

    throw AmiException(errorString.c_str());
}

void SteadystateProblem::writeErrorString(std::string *errorString,
                                          SteadyStateStatus status) const {
    /* write error message according to steady state status */
    switch (status) {
        case SteadyStateStatus::failed_too_long_simulation:
            (*errorString).append(": Simulation beyond t=1e100 without finding "
                                  " a steady state.");
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
    /* We need to check whether we still need to compute sensitivities.
       These boolean operation could be simplified, but here, clarity
       may more important than code reduction. */
    bool forwardSensisAlreadyComputed =
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        steady_state_status[1] == SteadyStateStatus::success &&
        model->getSteadyStateSensitivityMode() ==
        SteadyStateSensitivityMode::simulationFSA;
    bool needForwardSensisPosteq = !forwardSensisAlreadyComputed &&
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        solver->getSensitivityMethod() == SensitivityMethod::forward &&
        it > -1;
    bool needForwardSensisPreeq = !forwardSensisAlreadyComputed &&
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        it == -1;
    bool needForwardSensis = needForwardSensisPreeq || needForwardSensisPosteq;
    bool needForwardSensiByIntegration =
        solver->getSensitivityOrder() >= SensitivityOrder::first &&
        it == -1 && model->getSteadyStateSensitivityMode() ==
        SteadyStateSensitivityMode::simulationFSA;

    /* Check if we need to store sensis */
    switch (context) {
        case SteadyStateContext::newton:
            return needForwardSensis;
        case SteadyStateContext::storage:
            return needForwardSensis || forwardSensisAlreadyComputed;
        case SteadyStateContext::integration:
            return needForwardSensiByIntegration;
        default:
            throw AmiException("Requested invalid context in sensitivity "
                               "processing during steady state computation");
    }
}

realtype SteadystateProblem::getWrmsNorm(const AmiVector &x,
                                         const AmiVector &xdot,
                                         realtype atol,
                                         realtype rtol) {
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

bool SteadystateProblem::checkSteadyStateSuccess() const {
    /* Did one of the attempts yield s steady state? */
    if (std::any_of(steady_state_status.begin(), steady_state_status.end(),
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
            throw NewtonFailure(AMICI_NO_STEADY_STATE,
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
        const Solver *solver, Model *model, bool integrateForwardSensis) const
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
    if (integrateForwardSensis) {
        // need forward to compute sx0
        sim_solver->setSensitivityMethod(SensitivityMethod::forward);
    } else {
        sim_solver->setSensitivityMethod(SensitivityMethod::none);
        sim_solver->setSensitivityOrder(SensitivityOrder::none);
    }
    /* use x and sx as dummies for dx and sdx
     (they wont get touched in a CVodeSolver) */
    sim_solver->setup(model->t0(), model, x, x, sx, sx);

    return sim_solver;
}

void SteadystateProblem::getAdjointUpdates(Model &model,
                                           const ExpData &edata) {
    xB.reset();
    for (int it=0; it < model.nt(); it++) {
        if (std::isinf(model.getTimepoint(it))) {
            model.getAdjointStateObservableUpdate(
                slice(dJydx, it, model.nx_solver * model.nJ), it, x, edata);
            for (int ix = 0; ix < model.nxtrue_solver; ix++)
                xB[ix] += dJydx[ix + it * model.nx_solver];
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
