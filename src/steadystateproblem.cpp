#include "amici/steadystateproblem.h"
#include "amici/defines.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/solver_cvodes.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
#include "amici/newton_solver.h"
#include "amici/rdata.h"

#include <cmath>
#include <cstring>
#include <ctime>
#include <sundials/sundials_dense.h>
#include <memory>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_klu.h>

namespace amici {

void SteadystateProblem::workSteadyStateProblem(ReturnData *rdata,
                                               Solver *solver, Model *model,
                                               int it) {
    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param solver pointer to the AMICI solver object
     * @param model pointer to the AMICI model object
     * @param it integer with the index of the current time step
     * @param rdata pointer to the return data object
     */
    double run_time;
    clock_t starttime;

    /* First, try to do Newton steps */
    starttime = clock();

    auto newtonSolver = NewtonSolver::getSolver(t, x, solver->getLinearSolver(),
                                                model, rdata,
                                                solver->getNewtonMaxLinearSteps(),
                                                solver->getNewtonMaxSteps(),
                                                solver->getAbsoluteTolerance(),
                                                solver->getRelativeTolerance());

    auto newton_status = NewtonStatus::failed;
    try {
        applyNewtonsMethod(rdata, model, newtonSolver.get(), 1);
        newton_status = NewtonStatus::newt;
    } catch(NewtonFailure& ex) {
        try {
            /* Newton solver did not find a steady state, so try integration */
            getSteadystateSimulation(rdata, solver, model, it);
            newton_status = NewtonStatus::newt_sim;
        } catch(AmiException& ex) {
            /* may be integration failure from AmiSolve, so NewtonFailure
               won't do for all cases */
            try {
                applyNewtonsMethod(rdata, model, newtonSolver.get(), 2);
                newton_status = NewtonStatus::newt_sim_newt;
            } catch(NewtonFailure& ex) {
                throw amici::IntegrationFailure(ex.error_code,*t);
            }
        }
    } catch(...) {
        throw AmiException("Internal error in steady state problem");
    }
    run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;

    /* Compute steady state sensitvities */
    
    if (solver->getSensitivityOrder() >= SensitivityOrder::first &&
        (newton_status == NewtonStatus::newt ||
         newton_status == NewtonStatus::newt_sim_newt ||
         model->getSteadyStateSensitivityMode() != SteadyStateSensitivityMode::simulationFSA))
        // for newton_status == 2 the sensis were computed via FSA
        newtonSolver->computeNewtonSensis(sx);
    
    /* Get output of steady state solver, write it to x0 and reset time if necessary */
    writeNewtonOutput(rdata, model, newton_status, run_time, it);
}

realtype SteadystateProblem::getWrmsNorm(const AmiVector &x,const AmiVector &xdot) {
    N_VAbs(x.getNVector(), ewt.getNVector());
    N_VScale(rtol, ewt.getNVector(), ewt.getNVector());
    N_VAddConst(ewt.getNVector(), atol, ewt.getNVector());
    N_VInv(ewt.getNVector(), ewt.getNVector());
    return UNIT_ROUNDOFF * N_VWrmsNorm(xdot.getNVector(), ewt.getNVector());
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SteadystateProblem::applyNewtonsMethod(ReturnData *rdata,
                                           Model *model,
                                           NewtonSolver *newtonSolver,
                                           int newton_try) {
    int i_newtonstep = 0;
    int ix = 0;
    realtype wrms_tmp;
    double gamma = 1.0;
    bool compNewStep = TRUE;

    /* initialize output of linear solver for Newton step */
    delta.reset();

    /* Check, how fxdot is used exactly within AMICI... */
    model->fxdot(*t, x, &dx, &xdot);

    /* Check for relative error, but make sure not to divide by 0!
        Ensure positivity of the state */
    x_newton = *x;
    x_old = *x;
    xdot_old = xdot;
    
    //rdata->newton_numsteps[newton_try - 1] = 0.0;
    realtype wrms = getWrmsNorm(x_newton, xdot);
    bool converged = wrms < RCONST(1.0);
    while (!converged && i_newtonstep < newtonSolver->maxsteps) {

        /* If Newton steps are necessary, compute the inital search direction */
        if (compNewStep) {
            try{
                delta = xdot;
                newtonSolver->getStep(newton_try, i_newtonstep, &delta);
            } catch(NewtonFailure const& ex) {
                rdata->newton_numsteps[newton_try - 1] = static_cast<int>(getNaN());
                throw;
            } catch(std::exception const& ex) {
                rdata->newton_numsteps[newton_try - 1] = static_cast<int>(getNaN());
                throw NewtonFailure(AMICI_ERROR,"Newton method failed to compute new step!");
            }
        }
        
        /* Try a full, undamped Newton step */
        N_VLinearSum(1.0, x_old.getNVector(), gamma, delta.getNVector(), x->getNVector());
        
        /* Compute new xdot and residuals */
        model->fxdot(*t, x, &dx, &xdot);
        wrms_tmp = getWrmsNorm(x_newton, xdot);
        
        if (wrms_tmp < wrms) {
            /* If new residuals are smaller than old ones, update state */
            wrms = wrms_tmp;
            x_old = *x;
            xdot_old = xdot;
            /* New linear solve due to new state */
            compNewStep = TRUE;
            /* Check residuals vs tolerances */
            converged = wrms < RCONST(1.0);
            
            if (converged) {
                /* Ensure positivity of the found state */
                for (ix = 0; ix < model->nx; ix++) {
                    if ((*x)[ix] < 0.0) {
                        (*x)[ix] = 0.0;
                        converged = FALSE;
                    }
                }
            } else {
                /* increase dampening factor (superfluous, if converged) */
                gamma = fmin(1.0, 2.0 * gamma);
            }
        } else {
            /* Reduce dampening factor */
            gamma = gamma / 4.0;
            /* No new linear solve, only try new dampening */
            compNewStep = FALSE;
        }
        /* increase step counter */
        i_newtonstep++;
    }
    
    /* Set return values */
    rdata->newton_numsteps[newton_try-1] = i_newtonstep;
    if (!converged)
        throw NewtonFailure(AMICI_CONV_FAILURE,"applyNewtonsMethod");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SteadystateProblem::writeNewtonOutput(ReturnData *rdata,const Model *model,
                                         NewtonStatus newton_status,
                                         double run_time, int it)
{

    /* Get time for Newton solve */
    rdata->newton_time = run_time;
    rdata->newton_status = static_cast<int>(newton_status) ;
    if (newton_status == NewtonStatus::newt_sim) {
        rdata->t_steadystate = *t;
    }
    
    /* Steady state was found: set t to t0 if preeq, otherwise to inf */
    if (it == AMICI_PREEQUILIBRATE) {
        *t = model->t0();
    } else {
        *t = INFINITY;
    }
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SteadystateProblem::getSteadystateSimulation(ReturnData *rdata, Solver *solver,
                                                  Model *model, int it)
{
 
    std::unique_ptr<CVodeSolver> newtonSimSolver;

    /* Newton solver did not work, so try a simulation */
    if (it<1) {
        /* Preequilibration: Create a new CVode object for simulation */
        *t = model->t0();
        newtonSimSolver = createSteadystateSimSolver(solver, model, *t);
    } else {
        /* Carry on simulating from last point */
        *t = model->t(it-1);
    }
    
    /* Loop over steps and check for convergence */
    realtype wrms = getWrmsNorm(*x, xdot);
    bool converged = wrms < RCONST(1.0);
    
    int steps_newton = 0;
    while(!converged) {
        /* One step of ODE integration
         reason for tout specification:
         max with 1 ensures correct direction (any positive value would do)
         multiplication with 10 ensures nonzero difference and should ensure stable computation
         value is not important for AMICI_ONE_STEP mode, only direction w.r.t. current t
         */
        if (it<1)
            newtonSimSolver->solve(std::max(*t,1.0) * 10, x, &dx, t, AMICI_ONE_STEP);
        else
            solver->solve(std::max(*t,1.0) * 10, x, &dx, t, AMICI_ONE_STEP);

        /* Check for convergence */
        model->fxdot(*t, x, &dx, &xdot);
        wrms = getWrmsNorm(*x, xdot);
        converged = wrms < RCONST(1.0);
        /* increase counter, check for maxsteps */
        steps_newton++;
        if (steps_newton >= solver->getMaxSteps() && !converged)
            throw NewtonFailure(AMICI_TOO_MUCH_WORK,"getSteadystateSimulation");
    }
    if (it<1 && newtonSimSolver->getSensitivityOrder()>SensitivityOrder::none)
        newtonSimSolver->getSens(t, sx);
    else if (it>=0 && solver->getSensitivityOrder()>SensitivityOrder::none)
        solver->getSens(t, sx);
}

std::unique_ptr<CVodeSolver> SteadystateProblem::createSteadystateSimSolver(
        Solver *solver, Model *model, realtype tstart)
{

    
    /* Create new CVode object */
    
    auto newton_solver = std::unique_ptr<CVodeSolver>(new CVodeSolver());
    
    newton_solver->setLinearMultistepMethod(solver->getLinearMultistepMethod());
    newton_solver->setNonlinearSolverIteration(solver->getNonlinearSolverIteration());
    newton_solver->setAbsoluteTolerance(solver->getAbsoluteTolerance());
    newton_solver->setRelativeTolerance(solver->getRelativeTolerance());
    newton_solver->setMaxSteps(solver->getMaxSteps());
    newton_solver->setStabilityLimitFlag(solver->getStabilityLimitFlag());
    switch(solver->getLinearSolver()) {
        case LinearSolver::dense:
        case LinearSolver::KLU:
            newton_solver->setLinearSolver(solver->getLinearSolver());
            break;
        default:
            throw NewtonFailure(AMICI_NOT_IMPLEMENTED, "createSteadystateSimSolver");
    }
    newton_solver->setSensitivityOrder(solver->getSensitivityOrder());
    if (solver->getSensitivityMethod() != SensitivityMethod::none
        && model->getSteadyStateSensitivityMode() == SteadyStateSensitivityMode::simulationFSA)
        newton_solver->setSensitivityMethod(SensitivityMethod::forward); //need forward to compute sx0
    else
        newton_solver->setSensitivityMethod(SensitivityMethod::none);
    
    // use x and sx as dummies for dx and sdx (they wont get touched in a CVodeSolver)
    newton_solver->setup(x,x,sx,sx,model);
    
    return newton_solver;
}

} // namespace amici
