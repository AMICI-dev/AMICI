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
      dJydx(model.nJ * model.nx_solver * model.nt(), 0.0), numsteps(3, 0),
      numlinsteps(solver.getNewtonMaxSteps() * 2, 0) {}

void SteadystateProblem::workSteadyStateProblem(Solver *solver, Model *model,
                                                int it) {

    clock_t starttime;

    /* First, try to do Newton steps */
    starttime = clock();

    if (it == -1) {
        // solver was not run before, set up everything
        model->initialize(x, dx, sx, sdx,
                          solver->getSensitivityOrder() >=
                              SensitivityOrder::first);
        t = model->t0();
        solver->setup(t, model, x, dx, sx, sdx);
    } else {
        // solver was run before, extract current state from solver
        solver->writeSolution(&t, x, dx, sx);
    }

    auto newtonSolver = NewtonSolver::getSolver(&t, &x, *solver, model);

    newton_status = NewtonStatus::failed;
    try {
        applyNewtonsMethod(model, newtonSolver.get(), NewtonStatus::newt);
        newton_status = NewtonStatus::newt;
        std::copy_n(newtonSolver->getNumLinSteps().begin(),
                    solver->getNewtonMaxSteps(), numlinsteps.begin());
    } catch (NewtonFailure const &ex1) {
        std::copy_n(newtonSolver->getNumLinSteps().begin(),
                    solver->getNewtonMaxSteps(), numlinsteps.begin());
        try {
            /* Newton solver did not work, so try a simulation */
            if (it < 1) /* No previous time point computed, set t = t0 */
                t = model->t0();
            else /* Carry on simulating from last point */
                t = model->getTimepoint(it - 1);
            if (it < 0) {
                /* Preequilibration? -> Create a new CVode object for sim */
                auto newtonSimSolver =
                    createSteadystateSimSolver(solver, model);
                getSteadystateSimulation(newtonSimSolver.get(), model);
            } else {
                /* Solver was already created, use this one */
                getSteadystateSimulation(solver, model);
            }
            newton_status = NewtonStatus::newt_sim;
        } catch (AmiException const &ex2) {
            /* may be integration failure from AmiSolve, so NewtonFailure
               won't do for all cases */
            try {
                applyNewtonsMethod(model, newtonSolver.get(),
                                   NewtonStatus::newt_sim_newt);
                newton_status = NewtonStatus::newt_sim_newt;
                std::copy_n(
                    newtonSolver->getNumLinSteps().begin(),
                    solver->getNewtonMaxSteps(),
                    &numlinsteps.at(solver->getNewtonMaxLinearSteps() + 1));
            } catch (NewtonFailure const &ex3) {
                std::copy_n(
                    newtonSolver->getNumLinSteps().begin(),
                    solver->getNewtonMaxSteps(),
                    &numlinsteps.at(solver->getNewtonMaxLinearSteps() + 1));
                if (ex3.error_code == AMICI_TOO_MUCH_WORK)
                    throw AmiException("Steady state computation failed to "
                                       "converge within the allowed maximum "
                                       "number of iterations");
                throw;
            }
        }
    }
    cpu_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Compute steady state sensitvities */

    if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
        (newton_status == NewtonStatus::newt ||
         newton_status == NewtonStatus::newt_sim_newt ||
         model->getSteadyStateSensitivityMode() !=
             SteadyStateSensitivityMode::simulationFSA))
        // for newton_status == 2 the sensis were computed via FSA
        newtonSolver->computeNewtonSensis(sx);

    /* Get output of steady state solver, write it to x0 and reset time
     if necessary */
    storeSimulationState(model, solver->getSensitivityOrder() >=
                                    SensitivityOrder::first);
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
    model->fxdot(t, x, dx, xdot);
    wrms = getWrmsNorm(x, xdot, solver->getAbsoluteToleranceSteadyState(),
                       solver->getRelativeToleranceSteadyState());
    bool converged = wrms < RCONST(1.0);
    if (solver->getSensitivityOrder() > SensitivityOrder::none &&
        solver->getSensitivityMethod() == SensitivityMethod::forward) {
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

void SteadystateProblem::applyNewtonsMethod(Model *model,
                                            NewtonSolver *newtonSolver,
                                            NewtonStatus steadystate_try) {
    int i_newtonstep = 0;
    int ix = 0;
    double gamma = 1.0;
    bool compNewStep = TRUE;

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
                newtonSolver->getStep(steadystate_try == NewtonStatus::newt ? 1
                                                                            : 2,
                                      i_newtonstep, delta);
            } catch (NewtonFailure const &ex) {
                numsteps.at(steadystate_try == NewtonStatus::newt ? 0 : 2) =
                    i_newtonstep;
                throw;
            } catch (std::exception const &ex) {
                numsteps.at(steadystate_try == NewtonStatus::newt ? 0 : 2) =
                    i_newtonstep;
                throw AmiException("Newton solver failed to compute new step: "
                                   "%s",
                                   ex.what());
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
            compNewStep = TRUE;
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
              throw AmiException("Newton solver failed: a damping factor reached its lower bound");

            /* No new linear solve, only try new dampening */
            compNewStep = FALSE;
        }
        /* increase step counter */
        i_newtonstep++;
    }

    /* Set return values */
    numsteps.at(steadystate_try == NewtonStatus::newt ? 0 : 2) = i_newtonstep;
    if (!converged)
        throw NewtonFailure(AMICI_TOO_MUCH_WORK, "applyNewtonsMethod");
}

void SteadystateProblem::getSteadystateSimulation(Solver *solver,
                                                  Model *model)
{
    /* Loop over steps and check for convergence */
    bool converged = checkConvergence(solver, model);

    int steps_newton = 0;
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
        steps_newton++;
        if (steps_newton >= solver->getMaxSteps() && !converged) {
            numsteps.at(static_cast<int>(NewtonStatus::newt_sim) - 1) =
                steps_newton;
            throw NewtonFailure(AMICI_TOO_MUCH_WORK,
                                "exceeded maximum number of steps");
        }
    }
    numsteps.at(static_cast<int>(NewtonStatus::newt_sim) - 1) = steps_newton;
    if (solver->getSensitivityOrder() > SensitivityOrder::none)
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
