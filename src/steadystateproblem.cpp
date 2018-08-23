#include "amici/steadystateproblem.h"
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

    int newton_status;
    try {
        applyNewtonsMethod(rdata, model, newtonSolver.get(), 1);
        newton_status = 1;
    } catch(NewtonFailure& ex) {
        try {
            /* Newton solver did not find a steady state, so try integration */
            getSteadystateSimulation(rdata, solver, model, it);
            newton_status = 2;
        } catch(AmiException& ex) {// may be integration failure from AmiSolve, so NewtonFailure won't do for all cases
            try {
                applyNewtonsMethod(rdata, model, newtonSolver.get(), 2);
                newton_status = 3;
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
        solver->getSensitivityMethod() != SensitivityMethod::none) {
        // for newton_status == 2 the sensis were computed via FSA
        if (newton_status == 1 || newton_status == 3 || model->getSteadyStateSensitivityMode() != SteadyStateSensitivityMode::simulationFSA)
            newtonSolver->computeNewtonSensis(sx);
    
        if (it == AMICI_PREEQUILIBRATE) {
            for (int ip = 0; ip < model->nplist(); ip++) {
                for (int ix = 0; ix < model->nx; ix++) {
                    rdata->sx0[ip * model->nx + ix] = sx->at(ix,ip);
                }
            }
        }
    }
    
    
    

    /* Get output of steady state solver, write it to x0 and reset time if necessary */
    getNewtonOutput(rdata, model, newton_status, run_time, it);
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::applyNewtonsMethod(ReturnData *rdata,
                                           Model *model,
                                           NewtonSolver *newtonSolver,
                                           int newton_try) {
    /**
     * Runs the Newton solver iterations and checks for convergence to steady
     * state
     *
     * @param rdata pointer to the return data object
     * @param model pointer to the AMICI model object
     * @param newtonSolver pointer to the NewtonSolver object @type
     * NewtonSolver
     * @param newton_try integer start number of Newton solver (1 or 2)
     */
    int i_newtonstep = 0;
    int ix = 0;
    double res_tmp;
    double gamma = 1.0;
    bool compNewStep = TRUE;

    /* initialize output of linear solver for Newton step */
    delta.reset();

    /* Check, how fxdot is used exactly within AMICI... */
    model->fxdot(*t, x, &dx, &xdot);
    double res_abs = sqrt(N_VDotProd(xdot.getNVector(), xdot.getNVector()));

    /* Check for relative error, but make sure not to divide by 0!
        Ensure positivity of the state */
    x_newton = *x;
    N_VAbs(x_newton.getNVector(), x_newton.getNVector());
    N_VDiv(xdot.getNVector(), x_newton.getNVector(), rel_x_newton.getNVector());
    double res_rel = sqrt(N_VDotProd(rel_x_newton.getNVector(), rel_x_newton.getNVector()));
    x_old = *x;
    xdot_old = xdot;
    
    //rdata->newton_numsteps[newton_try - 1] = 0.0;
    bool converged = (res_abs < newtonSolver->atol || res_rel < newtonSolver->rtol);
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
        N_VDiv(xdot.getNVector(), x->getNVector(), rel_x_newton.getNVector());
        res_rel = sqrt(N_VDotProd(rel_x_newton.getNVector(), rel_x_newton.getNVector()));
        res_tmp = sqrt(N_VDotProd(xdot.getNVector(), xdot.getNVector()));
        
        if (res_tmp < res_abs) {
            /* If new residuals are smaller than old ones, update state */
            res_abs = res_tmp;
            x_old = *x;
            xdot_old = xdot;
            /* New linear solve due to new state */
            compNewStep = TRUE;
            /* Check residuals vs tolerances */
            converged = (res_abs < newtonSolver->atol) || (res_rel < newtonSolver->rtol);
            
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

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::getNewtonOutput(ReturnData *rdata,const Model *model,
                                         int newton_status,
                                         double run_time, int it) {
    /**
     * Stores output of workSteadyStateProblem in return data
     *
     * @param newton_status integer flag indicating when a steady state was
     * found
     * @param run_time double coputation time of the solver in milliseconds
     * @param rdata pointer to the return data instance
     * @param model pointer to the model instance
     * @param it current timepoint index, <0 indicates preequilibration
     */

    /* Get time for Newton solve */
    rdata->newton_time = run_time;
    rdata->newton_status = newton_status;
    
    /* Steady state was found: set t to t0 if preeq, otherwise to inf */
    if (it == AMICI_PREEQUILIBRATE) {
        *t = model->t0();

        /* Write steady state to output */
        rdata->x0 = x->getVector();
    } else {
        *t = INFINITY;
    }
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::getSteadystateSimulation(ReturnData *rdata, Solver *solver,
                                                  Model *model, int it) {
    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param solver pointer to the AMICI solver object
     * @param model pointer to the AMICI model object
     * @param rdata pointer to the return data object
     * @param it current timepoint index, <0 indicates preequilibration
     */
 
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
    double res_abs = INFINITY;
    double res_rel = INFINITY;
    
    int it_newton = 0;
    while(res_abs > solver->getAbsoluteTolerance() && res_rel > solver->getRelativeTolerance()) {
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

        model->fxdot(*t, x, &dx, &xdot);
        res_abs = sqrt(N_VDotProd(xdot.getNVector(), xdot.getNVector()));
        
        /* Ensure positivity and compute relative residual */
        x_newton = *x;
        N_VAbs(x_newton.getNVector(), x_newton.getNVector());
        for (int ix = 0; ix < model->nx; ix++)
            if (x_newton[ix] < solver->getAbsoluteTolerance())
                x_newton[ix] = solver->getAbsoluteTolerance();
        
        N_VDiv(xdot.getNVector(), x_newton.getNVector(), rel_x_newton.getNVector());
        res_rel = sqrt(N_VDotProd(rel_x_newton.getNVector(), rel_x_newton.getNVector()));
        
        /* Check for convergence */
        if (res_abs < solver->getAbsoluteTolerance() || res_rel < solver->getRelativeTolerance())
            break;
        /* increase counter, check for maxsteps */
        it_newton++;
        if (it_newton >= solver->getMaxSteps())
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
    /**
     * initialize CVodeSolver instance for preequilibration simulation
     *
     * @param solver pointer to the AMICI solver object
     * @param model pointer to the AMICI model object
     * @param tstart time point for starting Newton simulation
     * @return solver instance
     */
    
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
