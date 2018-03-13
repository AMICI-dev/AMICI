#include "amici/steadystateproblem.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
#include "amici/newton_solver.h"
#include "amici/rdata.h"

#include <cmath>
#include <cstring>
#include <ctime>
#include <sundials/sundials_dense.h>
#include <memory>

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
     * @param[in] solver pointer to the AMICI solver object @type Solver
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] it integer with the index of the current time step
     * @param[out] rdata pointer to the return data object @type ReturnData
     */
    double run_time;
    clock_t starttime;

    /* First, try to do Newton steps */
    starttime = clock();

    auto newtonSolver = std::unique_ptr<NewtonSolver>(NewtonSolver::getSolver(t, x, solver->getLinearSolver(),
                                                                              model, rdata,
                                                                              solver->getNewtonMaxLinearSteps(),
                                                                              solver->getNewtonMaxSteps(),
                                                                              solver->getAbsoluteTolerance(),
                                                                              solver->getRelativeTolerance()));
                                                      
    int newton_status;
    try {
        applyNewtonsMethod(rdata, model, newtonSolver.get(), 1);
        newton_status = 1;
    } catch(NewtonFailure& ex) {
        try {
            /* Newton solver did not find a steady state, so try integration */
            getNewtonSimulation(rdata, solver, model, it);
            newton_status = 2;
        } catch(AmiException& ex) {// may be integration failure from AmiSolve, so NewtonFailure won't do for all cases
            try {
                applyNewtonsMethod(rdata, model, newtonSolver.get(), 2);
                newton_status = 3;
            } catch(NewtonFailure& ex) {
                // TODO: more informative NewtonFailure to give more informative error code
                throw amici::IntegrationFailure(AMICI_CONV_FAILURE,*t);
            } catch(...) {
                throw AmiException("Internal error in steady state problem");
            }
        } catch(...) {
            throw AmiException("Internal error in steady state problem");
        }
    } catch(...) {
        throw AmiException("Internal error in steady state problem");
    }
    run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
    getNewtonOutput(rdata, model, newton_status, run_time, it);

    /* Compute steady state sensitvities */
    if (rdata->sensi_meth == AMICI_SENSI_FSA &&
        rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
        newtonSolver.get()->getSensis(it, sx);

    /* Reinitialize solver with preequilibrated state */
    if (it == AMICI_PREEQUILIBRATE) {
        solver->AMIReInit(*t, x, &dx);
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
            if (rdata->sensi_meth == AMICI_SENSI_FSA)
                solver->AMISensReInit(solver->getInternalSensitivityMethod(), sx, &sdx);
    }
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
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] newtonSolver pointer to the NewtonSolver object @type
     * NewtonSolver
     * @param[in] newton_try integer start number of Newton solver (1 or 2)
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
    for (ix = 0; ix < model->nx; ix++) {
        if (x_newton[ix] < newtonSolver->atol) {
            x_newton[ix] = newtonSolver->atol;
        }
    }
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
            } catch(std::exception ex) {
                rdata->newton_numsteps[newton_try - 1] = static_cast<int>(getNaN());
                throw NewtonFailure("Newton method failed to compute new step!");
            }
        }
        
        /* Try a full, undamped Newton step */
        N_VLinearSum(1.0, x_old.getNVector(), gamma, delta.getNVector(), x->getNVector());
        /* Ensure positivity of the state */
        for (ix = 0; ix < model->nx; ix++)
            if ((*x)[ix] < newtonSolver->atol)
                (*x)[ix] = newtonSolver->atol;
        
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
            /* increase dampening factor (superfluous, if converged) */
            gamma = fmin(1.0, 2.0 * gamma);
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
        throw NewtonFailure("Newton method failed to converge!");
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::getNewtonOutput(ReturnData *rdata,
                                         Model *model, int newton_status,
                                         double run_time, int it) {
    /**
     * Stores output of workSteadyStateProblem in return data
     *
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] newton_status integer flag indicating when a steady state was
     * found
     * @param[in] run_time double coputation time of the solver in milliseconds
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[in] it current timepoint index, <0 indicates preequilibration @type int
     */

    /* Get time for Newton solve */
    rdata->newton_time = run_time;
    rdata->newton_status = newton_status;
    
    /* Steady state was found: set t to t0 if preeq, otherwise to inf */
    if (it == AMICI_PREEQUILIBRATE) {
        *t = rdata->ts[0];

        /* Write steady state to output */
        for (int ix = 0; ix < model->nx; ix++) {
            rdata->x0[ix] = (*x)[ix];
        }
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

void SteadystateProblem::getNewtonSimulation(ReturnData *rdata, Solver *solver,
                                            Model *model, int it) {
    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param[in] solver pointer to the AMICI solver object @type Solver
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[in] it current timepoint index, <0 indicates preequilibration @type int
     */
 
    realtype tstart;
    
    /* Newton solver did not work, so try a simulation: reinitialize solver */
    if (it<1)
        tstart = model->t0();
    else
        tstart = rdata->ts[it-1];
    *t = tstart;
    
    model->fx0(x);
    solver->AMIReInit(*t, x, &dx);
    
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
        solver->AMISolve(std::max(*t,1.0) * 10, x, &dx, t,
                                  AMICI_ONE_STEP);
        model->fxdot(*t, x, &dx, &xdot);
        res_abs = sqrt(N_VDotProd(xdot.getNVector(), xdot.getNVector()));
        
        /* Ensure positivity and compute relative residual */
        x_newton = *x;
        N_VAbs(x_newton.getNVector(), x_newton.getNVector());
        for (int ix = 0; ix < model->nx; ix++) {
            if (x_newton[ix] < solver->getAbsoluteTolerance()) {
                x_newton[ix] = solver->getAbsoluteTolerance();
            }
        }
        N_VDiv(xdot.getNVector(), x_newton.getNVector(), rel_x_newton.getNVector());
        res_rel = sqrt(N_VDotProd(rel_x_newton.getNVector(), rel_x_newton.getNVector()));
        
        /* Check for convergence */
        if (res_abs < solver->getAbsoluteTolerance() || res_rel < solver->getRelativeTolerance())
            break;
        /* increase counter, check for maxsteps */
        it_newton++;
        if (it_newton >= solver->getMaxSteps())
            throw NewtonFailure("Simulation based steady state failed to converge");
            
        
    }
}

} // namespace amici
