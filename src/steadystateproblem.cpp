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

SteadystateProblem::SteadystateProblem() {}

int SteadystateProblem::workSteadyStateProblem(UserData *udata, TempData *tdata,
                                               ReturnData *rdata, Solver *solver,
                                               Model *model, int it) {
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
     * @return stats integer flag indicating success of the method
     */

    int status = (int)*rdata->status;
    double run_time;
    clock_t starttime;

    /* First, try to do Newton steps */
    starttime = clock();

    NewtonSolver *newtonSolver = NewtonSolver::getSolver(
        udata->linsol, model, rdata, udata, tdata, &status);
    
    status = applyNewtonsMethod(udata, rdata, tdata, model, newtonSolver, 1);

    if (status == AMICI_SUCCESS) {
        /* if the Newton solver found a steady state */
        run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
        getNewtonOutput(tdata, rdata, model, 1, run_time, it);
    } else {
        /* Newton solver did not find a steady state, so try integration */
        status = getNewtonSimulation(udata, tdata, rdata, solver, model);

        if (status == AMICI_SUCCESS) {
            /* if simulation found a steady state */
            run_time = (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
            getNewtonOutput(tdata, rdata, model, 2, run_time, it);
        } else {
            status =
                applyNewtonsMethod(udata, rdata, tdata, model, newtonSolver, 2);

            if (status == AMICI_SUCCESS) {
                /* If the second Newton solver found a steady state */
                run_time =
                    (double)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
                getNewtonOutput(tdata, rdata, model, 3, run_time, it);
            }
        }
    }
    
    /* Compute steady state sensitvities */
    if (rdata->sensi_meth == AMICI_SENSI_FSA && rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
        if (status == AMICI_SUCCESS)
            status = newtonSolver->getSensis(it);
    
    /* Reinitialize solver with preequilibrated state */
    if (it == AMICI_PREEQUILIBRATE) {

        status = solver->AMIReInit(tdata->t, tdata->x, tdata->dx);
        if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST)
            if (rdata->sensi_meth == AMICI_SENSI_FSA)
                status = solver->AMISensReInit(udata->ism, tdata->sx, tdata->sdx);
        if (status != AMICI_SUCCESS)
            status = AMICI_ERROR_PREEQUILIBRATION;
    }

    delete newtonSolver;
    return status;
}

/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

int SteadystateProblem::applyNewtonsMethod(UserData *udata, ReturnData *rdata,
                                           TempData *tdata, Model *model,
                                           NewtonSolver *newtonSolver, int newton_try) {
    /**
     * Runs the Newton solver iterations and checks for convergence to steady state
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @param[out] tdata pointer to the temporary data object @type TempData
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] newtonSolver pointer to the NewtonSolver object @type NewtonSolver
     * @param[in] newton_try integer start number of Newton solver (1 or 2)
     * @return stats integer flag indicating success of the method
     */

    int status = AMICI_ERROR_NEWTONSOLVER;
    int i_newtonstep = 0;
    int ix = 0;
    double res_abs;
    double res_rel;
    double res_tmp;
    double gamma = 1.0;
    realtype *x_tmp;

    N_Vector delta = N_VNew_Serial(model->nx);
    N_Vector rel_x_newton = N_VNew_Serial(model->nx);
    N_Vector x_newton = N_VNew_Serial(model->nx);

    /* initialize output von linear solver for Newton step */
    N_VConst(0.0, delta);

    /* Check, how fxdot is used exactly within AMICI... */
    model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);
    res_abs = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));

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
    res_rel = sqrt(N_VDotProd(rel_x_newton, rel_x_newton));

    if (res_abs >= udata->atol && res_rel >= udata->rtol) {

        /* If Newton steps are necessary, compute the inital search direction */
        status = newtonSolver->getStep(newton_try, i_newtonstep, delta);

        if (status == AMICI_SUCCESS) {
            /* The linear solver was successful, now the Newton solver needs to
             * be */
            status = AMICI_ERROR_NEWTONSOLVER;

            /* Copy the current state to the old one, make up a new vector for
             * JDiag */
            N_VScale(1.0, tdata->x, tdata->x_old);
            N_VScale(1.0, tdata->xdot, tdata->xdot_old);

            /* Newton iterations */
            for (i_newtonstep = 0; i_newtonstep < udata->newton_maxsteps;
                 i_newtonstep++) {

                /* Try a full, undamped Newton step */
                N_VLinearSum(1.0, tdata->x_old, gamma, delta, tdata->x);

                /* Ensure positivity of the state */
                x_tmp = N_VGetArrayPointer(tdata->x);
                for (ix = 0; ix < model->nx; ix++) {
                    if (x_tmp[ix] < 0.0) {
                        x_tmp[ix] = 0.0;
                    }
                }

                /* Compute new xdot */
                model->fxdot(tdata->t, tdata->x, tdata->dx, tdata->xdot, tdata);

                /* Check if new residuals are smaller than old ones */
                res_tmp = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));

                if (res_tmp < res_abs) {
                    /* update state */
                    res_abs = res_tmp;
                    N_VScale(1.0, tdata->x, tdata->x_old);
                    N_VScale(1.0, tdata->xdot, tdata->xdot_old);

                    /* Check residuals vs tolerances */
                    if (res_abs < udata->atol) {
                        /* Return number of Newton steps */
                        rdata->newton_numsteps[newton_try - 1] =
                            i_newtonstep + 1;
                        status = AMICI_SUCCESS;
                        break;
                    }

                    if (status != AMICI_SUCCESS) {
                        /* increase dampening factor */
                        gamma = fmax(1.0, 2.0 * gamma);

                        /* Do another Newton step */
                        status = newtonSolver->getStep(newton_try, i_newtonstep,
                                                       delta);
                        if (status == AMICI_SUCCESS) {
                            /* Newton step was successful, now Newtons method
                             * still needs to be */
                            status = AMICI_ERROR_NEWTONSOLVER;
                        } else {
                            /* Linear solver errored, go to clean up and return
                             * part */
                            rdata->newton_numsteps[newton_try - 1] =
                                amiGetNaN();
                            break;
                        }
                    }
                } else {
                    /* Reduce dampening factor */
                    gamma = gamma / 4.0;
                }
            }

            /* Set return values */
            rdata->newton_numsteps[newton_try - 1] = i_newtonstep;
        } else {
            rdata->newton_numsteps[newton_try - 1] = amiGetNaN();
        }

    } else {
        /* No Newton steps were necessary */
        status = AMICI_SUCCESS;

        /* Set return values */
        rdata->newton_numsteps[newton_try - 1] = 0.0;
    }

    /* Clean up worksapce */
    N_VDestroy_Serial(delta);
    N_VDestroy_Serial(rel_x_newton);
    N_VDestroy_Serial(x_newton);

    return (status);
}

/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

void SteadystateProblem::getNewtonOutput(TempData *tdata, ReturnData *rdata,
                                        Model *model, int newton_status,
                                        double run_time, int it) {
    /**
     * Stores output of workSteadyStateProblem in return data
     *
     * @param[in] tdata pointer to the temporary data object @type UserData
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] newton_status integer flag indicating when a steady state was found
     * @param[in] run_time double coputation time of the solver in milliseconds
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @return stats integer flag indicating success of the method
     */

    /* Get time for Newton solve */
    rdata->newton_time[0] = run_time;

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

    /* Write Newton flag */
    *rdata->newton_status = (double)newton_status;
}

/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

int SteadystateProblem::getNewtonSimulation(UserData *udata, TempData *tdata,
                                            ReturnData *rdata, Solver *solver,
                                            Model *model) {
    /**
     * Forward simulation is launched, if Newton solver fails in first try
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] solver pointer to the AMICI solver object @type Solver
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[out] tdata pointer to the temporary data object @type TempData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @return stats integer flag indicating success of the method
     */
 
    double res_abs;
    double res_rel;
    double sim_time;
    int status = (int)*rdata->status;
    realtype *x_tmp;
    N_Vector rel_x_newton = N_VNew_Serial(model->nx);
    N_Vector x_newton = N_VNew_Serial(model->nx);

    /* Newton solver did not work, so try a simulation */
    if (tdata->t >= 1e6) {
        sim_time = 10.0 * (tdata->t);
    } else {
        sim_time = 1e6;
    }
    status = solver->AMISolve(RCONST(sim_time), tdata->x, tdata->dx,
                              &(tdata->t), AMICI_NORMAL);

    if (status == AMICI_SUCCESS) {
        /* Check residuals */
        res_abs = sqrt(N_VDotProd(tdata->xdot, tdata->xdot));

        /* Ensure positivity for relative residual */
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

        /* residuals are small? */
        if (res_abs < udata->atol || res_rel < udata->rtol) {
            return (AMICI_SUCCESS);
        } else {
            status = AMICI_ERROR_SIM2STEADYSTATE;
        }
    }

    N_VDestroy_Serial(rel_x_newton);
    N_VDestroy_Serial(x_newton);
    return status;
}

/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------- */

int SteadystateProblem::linsolveSPBCG(UserData *udata, ReturnData *rdata,
                                      TempData *tdata, Model *model, int ntry,
                                      int nnewt, N_Vector ns_delta) {
    /**
     * Iterative linear solver created from SPILS BiCG-Stab.
     * Solves the linear system within each Newton step if iterative solver is chosen.
     *
     * @param[in] udata pointer to the user data object @type UserData
     * @param[in] model pointer to the AMICI model object @type Model
     * @param[in] ntry integer newton_try integer start number of Newton solver (1 or 2)
     * @param[in] nnewt integer number of current Newton step
     * @param[in] ns_delta ???
     * @param[out] tdata pointer to the temporary data object @type TempData
     * @param[out] rdata pointer to the return data object @type ReturnData
     * @return stats integer flag indicating success of the method
     */

    int status = AMICI_ERROR_NEWTON_LINSOLVER;

    double rho;
    double rho1;
    double alpha;
    double beta;
    double omega;
    double res;

    N_Vector ns_p = N_VNew_Serial(model->nx);
    N_Vector ns_h = N_VNew_Serial(model->nx);
    N_Vector ns_t = N_VNew_Serial(model->nx);
    N_Vector ns_s = N_VNew_Serial(model->nx);
    N_Vector ns_r = N_VNew_Serial(model->nx);
    N_Vector ns_rt = N_VNew_Serial(model->nx);
    N_Vector ns_v = N_VNew_Serial(model->nx);
    N_Vector ns_Jv = N_VNew_Serial(model->nx);
    N_Vector ns_tmp = N_VNew_Serial(model->nx);
    N_Vector ns_Jdiag = N_VNew_Serial(model->nx);

    N_VScale(-1.0, tdata->xdot, tdata->xdot);

    // Get the diagonal of the Jacobian for preconditioning
    model->fJDiag(tdata->t, ns_Jdiag, tdata->x, tdata);

    // Ensure positivity of entries in ns_Jdiag
    N_VConst(1.0, ns_p);
    N_VAbs(ns_Jdiag, ns_tmp);
    N_VCompare(1e-15, ns_tmp, ns_tmp);
    N_VLinearSum(-1.0, ns_tmp, 1.0, ns_p, ns_tmp);
    N_VLinearSum(1.0, ns_Jdiag, 1.0, ns_tmp, ns_Jdiag);

    // Initialize for linear solve
    N_VConst(0.0, ns_p);
    N_VConst(0.0, ns_v);
    N_VConst(0.0, ns_delta);
    N_VConst(0.0, ns_tmp);
    rho = 1.0;
    omega = 1.0;
    alpha = 1.0;

    // can be set to 0 at the moment
    model->fJv(ns_delta, ns_Jv, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);

    // ns_r = xdot - ns_Jv;
    N_VLinearSum(-1.0, ns_Jv, 1.0, tdata->xdot, ns_r);
    N_VDiv(ns_r, ns_Jdiag, ns_r);
    res = sqrt(N_VDotProd(ns_r, ns_r));
    N_VScale(1.0, ns_r, ns_rt);

    for (int i_linstep = 0; i_linstep < udata->newton_maxlinsteps;
         i_linstep++) {
        // Compute factors
        rho1 = rho;
        rho = N_VDotProd(ns_rt, ns_r);
        beta = rho * alpha / (rho1 * omega);

        // ns_p = ns_r + beta * (ns_p - omega * ns_v);
        N_VLinearSum(1.0, ns_p, -omega, ns_v, ns_p);
        N_VLinearSum(1.0, ns_r, beta, ns_p, ns_p);

        // ns_v = J * ns_p
        model->fJv(ns_p, ns_v, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);
        N_VDiv(ns_v, ns_Jdiag, ns_v);

        // Compute factor
        alpha = rho / N_VDotProd(ns_rt, ns_v);

        // ns_h = ns_delta + alpha * ns_p;
        N_VLinearSum(1.0, ns_delta, alpha, ns_p, ns_h);
        // ns_s = ns_r - alpha * ns_v;
        N_VLinearSum(1.0, ns_r, -alpha, ns_v, ns_s);

        // ns_t = J * ns_s
        model->fJv(ns_s, ns_t, tdata->t, tdata->x, tdata->xdot, tdata, ns_tmp);
        N_VDiv(ns_t, ns_Jdiag, ns_t);

        // Compute factor
        omega = N_VDotProd(ns_t, ns_s) / N_VDotProd(ns_t, ns_t);

        // ns_delta = ns_h + omega * ns_s;
        N_VLinearSum(1.0, ns_h, omega, ns_s, ns_delta);
        // ns_r = ns_s - omega * ns_t;
        N_VLinearSum(1.0, ns_s, -omega, ns_t, ns_r);

        // Compute the (unscaled) residual
        N_VProd(ns_r, ns_Jdiag, ns_r);
        res = sqrt(N_VDotProd(ns_r, ns_r));

        // Test convergence
        if (res < udata->atol) {
            // Write number of steps needed
            rdata->newton_numlinsteps[(ntry - 1) * udata->newton_maxsteps +
                                      nnewt] = i_linstep + 1;

            // Return success
            N_VScale(-1.0, tdata->xdot, tdata->xdot);
            status = AMICI_SUCCESS;
            break;
        }

        // Scale back
        N_VDiv(ns_r, ns_Jdiag, ns_r);
    }

    // Clean up workspace
    N_VDestroy_Serial(ns_p);
    N_VDestroy_Serial(ns_h);
    N_VDestroy_Serial(ns_t);
    N_VDestroy_Serial(ns_s);
    N_VDestroy_Serial(ns_r);
    N_VDestroy_Serial(ns_rt);
    N_VDestroy_Serial(ns_v);
    N_VDestroy_Serial(ns_Jv);
    N_VDestroy_Serial(ns_tmp);
    N_VDestroy_Serial(ns_Jdiag);

    N_VScale(-1.0, tdata->xdot, tdata->xdot);

    // Return
    return (status);
}
