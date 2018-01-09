#include "../include/steadystateproblem.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include "include/edata.h"
#include "include/forwardproblem.h"
#include "include/newton_solver.h"
#include "include/rdata.h"
#include "include/tdata.h"
#include "include/udata.h"
#include <cstring>
#include <ctime>
#include <sundials/sundials_dense.h>
#include <memory>

namespace amici {

void SteadystateProblem::workSteadyStateProblem(const UserData *udata,
                                               TempData *tdata,
                                               ReturnData *rdata,
                                               Solver *solver, Model *model,
                                               int it) {
    /**
     * Tries to determine the steady state of the ODE system by a Newton
     * solver, uses forward intergration, if the Newton solver fails,
     * restarts Newton solver, if integration fails.
     * Computes steady state sensitivities
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] solver pointer to the AMICI solver object @type Solver
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] it integer with the index of the current time step
     * @param[out] tdata pointer to the temporary data object @type TempData
     * @param[out] rdata pointer to the return data object @type ReturnData
     */
    double run_time;
    clock_t starttime;

    /* First, try to do Newton steps */
    starttime = clock();

    auto newtonSolver = std::unique_ptr<NewtonSolver>(NewtonSolver::getSolver(udata->linsol, model, rdata, udata, tdata));
                                                      
    int newton_status;
    try {
        applyNewtonsMethod(udata, rdata, tdata, model, newtonSolver.get(), 1);
        newton_status = 1;
    } catch(NewtonFailure& ex) {
        try {
            /* Newton solver did not find a steady state, so try integration */
            getNewtonSimulation(udata, tdata, rdata, solver, model, it);
            newton_status = 2;
        } catch(AmiException& ex) {// may be integration failure from AmiSolve, so NewtonFailure won't do for all cases
            try {
                applyNewtonsMethod(udata, rdata, tdata, model, newtonSolver.get(), 2);
                newton_status = 3;
            } catch(NewtonFailure& ex) {
                // TODO: more informative NewtonFailure to give more informative error code
                throw amici::IntegrationFailure(AMICI_CONV_FAILURE,tdata->t);
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
    getNewtonOutput(tdata, rdata, model, newton_status, run_time, it);

    /* Compute steady state sensitvities */
    if (rdata->sensi_meth == AMICI_SENSI_FSA &&
        rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
        newtonSolver->getSensis(it);

    /* Reinitialize solver with preequilibrated state */
    if (it == AMICI_PREEQUILIBRATE) {
        solver->AMIReInit(tdata->t, tdata->x, tdata->dx);
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
            if (rdata->sensi_meth == AMICI_SENSI_FSA)
                solver->AMISensReInit(udata->ism, tdata->sx, tdata->sdx);
    }
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::applyNewtonsMethod(const UserData *udata,
                                           ReturnData *rdata, TempData *tdata,
                                           Model *model,
                                           NewtonSolver *newtonSolver,
                                           int newton_try) {
    /**
     * Runs the Newton solver iterations and checks for convergence to steady
     * state
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[out] tdata pointer to the temporary data object @type TempData
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
    
    realtype *x_tmp;

    /* initialize output von linear solver for Newton step */
    N_VConst(0.0, delta);

    /* Check, how fxdot is used exactly within AMICI... */
    model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
    double res_abs = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));

    /* Check for relative error, but make sure not to divide by 0!
        Ensure positivity of the state */
    N_VScale(1.0, tdata->x, x_newton);
    N_VAbs(x_newton, x_newton);
    x_tmp = N_VGetArrayPointer(x_newton);
    for (ix = 0; ix < model->nx; ix++) {
        if (x_tmp[ix] < udata->atol) {
            x_tmp[ix] = udata->atol;
        }
    }
    N_VDiv(tdata->xdot, x_newton, rel_x_newton);
    double res_rel = sqrt(N_VDotProd(rel_x_newton, rel_x_newton));
    N_VScale(1.0, tdata->x, tdata->x_old);
    N_VScale(1.0, tdata->xdot, tdata->xdot_old);
    
    //rdata->newton_numsteps[newton_try - 1] = 0.0;
    bool converged = (res_abs < udata->atol || res_rel < udata->rtol);
    while (!converged && i_newtonstep < udata->newton_maxsteps) {

        /* If Newton steps are necessary, compute the inital search direction */
        if (compNewStep) {
            try{
                newtonSolver->getStep(newton_try, i_newtonstep, delta);
            } catch(...) {
                rdata->newton_numsteps[newton_try - 1] = amiGetNaN();
                throw NewtonFailure("Newton method failed to compute new step!");
            }
        }
        
        /* Try a full, undamped Newton step */
        N_VLinearSum(1.0, tdata->x_old, gamma, delta, tdata->x);
        /* Ensure positivity of the state */
        x_tmp = N_VGetArrayPointer(tdata->x);
        for (ix = 0; ix < model->nx; ix++)
            if (x_tmp[ix] < udata->atol)
                x_tmp[ix] = udata->atol;
        
        /* Compute new xdot and residuals */
        model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
        N_VDiv(tdata->xdot, tdata->x, rel_x_newton);
        res_rel = sqrt(N_VDotProd(rel_x_newton, rel_x_newton));
        res_tmp = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));
        
        if (res_tmp < res_abs) {
            /* If new residuals are smaller than old ones, update state */
            res_abs = res_tmp;
            N_VScale(1.0, tdata->x, tdata->x_old);
            N_VScale(1.0, tdata->xdot, tdata->xdot_old);
            /* New linear solve due to new state */
            compNewStep = TRUE;
            /* Check residuals vs tolerances */
            converged = (res_abs < udata->atol) || (res_rel < udata->rtol);
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
    rdata->newton_numsteps[newton_try-1] = (double) i_newtonstep;
    if (!converged)
        throw NewtonFailure("Newton method failed to converge!");
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::getNewtonOutput(TempData *tdata, ReturnData *rdata,
                                         Model *model, int newton_status,
                                         double run_time, int it) {
    /**
     * Stores output of workSteadyStateProblem in return data
     *
     * @param[in] tdata pointer to the temporary data object @type UserData
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] newton_status integer flag indicating when a steady state was
     * found
     * @param[in] run_time double coputation time of the solver in milliseconds
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[in] it current timepoint index, <0 indicates preequilibration @type int
     */

    /* Get time for Newton solve */
    rdata->newton_time[0] = run_time;
    rdata->newton_status[0] = newton_status;
    
    /* Steady state was found: set t to t0 if preeq, otherwise to inf */
    if (it == AMICI_PREEQUILIBRATE) {
        tdata->t = rdata->ts[0];

        /* Write steady state to output */
        realtype *x_tmp = NV_DATA_S(tdata->x);
        for (int ix = 0; ix < model->nx; ix++) {
            rdata->x0[ix] = x_tmp[ix];
        }
    } else {
        tdata->t = INFINITY;
    }
}

/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------
 */

void SteadystateProblem::getNewtonSimulation(const UserData *udata, TempData *tdata,
                                            ReturnData *rdata, Solver *solver,
                                            Model *model, int it) {
    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] solver pointer to the AMICI solver object @type Solver
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[out] tdata pointer to the temporary data object @type TempData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[in] it current timepoint index, <0 indicates preequilibration @type int
     */
 
    realtype tstart;
    
    /* Newton solver did not work, so try a simulation: reinitialize solver */
    if (it<1)
        tstart = udata->tstart;
    else
        tstart = udata->ts[it-1];
    tdata->t = tstart;
    
    model->fx0(tdata->x, tdata);
    solver->AMIReInit(udata->tstart, tdata->x, tdata->dx);
    
    /* Loop over steps and check for convergence */
    double res_abs = INFINITY;
    double res_rel = INFINITY;
    realtype *x_tmp;
    
    int it_newton = 0;
    while(res_abs > udata->atol && res_rel > udata->rtol) {
        /* One step of ODE integration */
        solver->AMISolve(1e12, tdata->x, tdata->dx, &(tdata->t),
                                  AMICI_ONE_STEP);
        model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
        res_abs = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));
        
        /* Ensure positivity and compute relative residual */
        N_VScale(1.0, tdata->x, x_newton);
        N_VAbs(x_newton, x_newton);
        x_tmp = N_VGetArrayPointer(x_newton);
        for (int ix = 0; ix < model->nx; ix++) {
            if (x_tmp[ix] < udata->atol) {
                x_tmp[ix] = udata->atol;
            }
        }
        N_VDiv(tdata->xdot, x_newton, rel_x_newton);
        res_rel = sqrt(N_VDotProd(rel_x_newton, rel_x_newton));
        
        /* Check for convergence */
        if (res_abs < udata->atol || res_rel < udata->rtol) 
            break;
        /* increase counter, check for maxsteps */
        it_newton++;
        if (it_newton >= udata->maxsteps)
            throw NewtonFailure("Simulation based steady state failed to converge");
            
        
    }
}

} // namespace amici