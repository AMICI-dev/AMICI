#include "amici/solver_cvodes.h"

#include "amici/exception.h"
#include "amici/model_ode.h"
#include "amici/sundials_linsol_wrapper.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_impl.h>

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#include <sstream>

#define ZERO SUN_RCONST(0.0)
#define ONE SUN_RCONST(1.0)
#define FOUR SUN_RCONST(4.0)

namespace amici {

// Ensure AMICI options and return codes are in sync with SUNDIALS
static_assert((int)InternalSensitivityMethod::simultaneous == CV_SIMULTANEOUS);
static_assert((int)InternalSensitivityMethod::staggered == CV_STAGGERED);
static_assert((int)InternalSensitivityMethod::staggered1 == CV_STAGGERED1);

static_assert((int)InterpolationType::hermite == CV_HERMITE);
static_assert((int)InterpolationType::polynomial == CV_POLYNOMIAL);

static_assert((int)LinearMultistepMethod::adams == CV_ADAMS);
static_assert((int)LinearMultistepMethod::BDF == CV_BDF);

#define STATIC_ASSERT_EQUAL(amici_constant, cv_constant)                       \
    static_assert(                                                             \
        amici_constant == cv_constant, #amici_constant " != " #cv_constant     \
    )

STATIC_ASSERT_EQUAL(amici::AMICI_SUCCESS, CV_SUCCESS);
STATIC_ASSERT_EQUAL(amici::AMICI_ROOT_RETURN, CV_ROOT_RETURN);
STATIC_ASSERT_EQUAL(amici::AMICI_DATA_RETURN, CV_TSTOP_RETURN);
STATIC_ASSERT_EQUAL(amici::AMICI_ILL_INPUT, CV_ILL_INPUT);
STATIC_ASSERT_EQUAL(amici::AMICI_NORMAL, CV_NORMAL);
STATIC_ASSERT_EQUAL(amici::AMICI_ONE_STEP, CV_ONE_STEP);
STATIC_ASSERT_EQUAL(amici::AMICI_TOO_MUCH_ACC, CV_TOO_MUCH_ACC);
STATIC_ASSERT_EQUAL(amici::AMICI_TOO_MUCH_WORK, CV_TOO_MUCH_WORK);
STATIC_ASSERT_EQUAL(amici::AMICI_ERR_FAILURE, CV_ERR_FAILURE);
STATIC_ASSERT_EQUAL(amici::AMICI_CONV_FAILURE, CV_CONV_FAILURE);
STATIC_ASSERT_EQUAL(amici::AMICI_LSETUP_FAIL, CV_LSETUP_FAIL);
STATIC_ASSERT_EQUAL(amici::AMICI_CONSTR_FAIL, CV_CONSTR_FAIL);
STATIC_ASSERT_EQUAL(amici::AMICI_RHSFUNC_FAIL, CV_RHSFUNC_FAIL);
STATIC_ASSERT_EQUAL(amici::AMICI_FIRST_RHSFUNC_ERR, CV_FIRST_RHSFUNC_ERR);
STATIC_ASSERT_EQUAL(amici::AMICI_FIRST_QRHSFUNC_ERR, CV_FIRST_QRHSFUNC_ERR);
STATIC_ASSERT_EQUAL(amici::AMICI_BAD_T, CV_BAD_T);
STATIC_ASSERT_EQUAL(amici::AMICI_BAD_DKY, CV_BAD_DKY);
STATIC_ASSERT_EQUAL(amici::AMICI_SRHSFUNC_FAIL, CV_SRHSFUNC_FAIL);
STATIC_ASSERT_EQUAL(amici::AMICI_FIRST_SRHSFUNC_ERR, CV_FIRST_SRHSFUNC_ERR);
STATIC_ASSERT_EQUAL(amici::AMICI_REPTD_SRHSFUNC_ERR, CV_REPTD_SRHSFUNC_ERR);
STATIC_ASSERT_EQUAL(amici::AMICI_UNREC_SRHSFUNC_ERR, CV_UNREC_SRHSFUNC_ERR);
STATIC_ASSERT_EQUAL(amici::AMICI_RTFUNC_FAIL, CV_RTFUNC_FAIL);

/*
 * The following static members are callback function to CVODES.
 * Their signatures must not be changes.
 */
static int fxdot(realtype t, N_Vector x, N_Vector xdot, void* user_data);

static int fJSparse(
    realtype t, N_Vector x, N_Vector xdot, SUNMatrix J, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
);

static int
fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J, void* user_data,
   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int
fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

static int fJSparseB(
    realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B
);

static int fJBand(
    realtype t, N_Vector x, N_Vector xdot, SUNMatrix J, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
);

static int fJBandB(
    realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B
);

static int
fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot,
    void* user_data, N_Vector tmp);

static int fJvB(
    N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB,
    N_Vector xBdot, void* user_data, N_Vector tmpB
);

static int froot(realtype t, N_Vector x, realtype* root, void* user_data);

static int
fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void* user_data);

static int
fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void* user_data);

static int fxBdot_ss(realtype t, N_Vector xB, N_Vector xBdot, void* user_data);

static int fqBdot_ss(realtype t, N_Vector xB, N_Vector qBdot, void* user_data);

static int fJSparseB_ss(
    realtype t, N_Vector x, N_Vector xBdot, SUNMatrix JB, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
);

static int fsxdot(
    int Ns, realtype t, N_Vector x, N_Vector xdot, int ip, N_Vector sx,
    N_Vector sxdot, void* user_data, N_Vector tmp1, N_Vector tmp2
);

/* Function implementations */

void CVodeSolver::init(
    realtype const t0, AmiVector const& x0, AmiVector const& /*dx0*/
) const {
    solver_was_called_F_ = false;
    force_reinit_postprocess_F_ = false;
    t_ = t0;
    x_ = x0;
    int status;
    if (get_init_done()) {
        status = CVodeReInit(solver_memory_.get(), t0, x_.get_nvector());
    } else {
        status = CVodeInit(solver_memory_.get(), fxdot, t0, x_.get_nvector());
        set_init_done();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeInit");
}

void CVodeSolver::init_steady_state(
    realtype const /*t0*/, AmiVector const& /*x0*/, AmiVector const& /*dx0*/
) const {
    // We need to set the steadystate rhs function. Sundials doesn't have this
    // in its public API, so we have to change it in the solver memory,
    // as re-calling init would unset solver settings.
    auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
    cv_mem->cv_f = fxBdot_ss;

    // Since SUNDIALS v5.8.0, we also need to update the NlsRhs function,
    // otherwise the old value of `cv_mem->cv_f` would still be used there and
    // lead to incorrect simulation results.
    CVodeSetNlsRhsFn(solver_memory_.get(), fxBdot_ss);
}

void CVodeSolver::sens_init_1(
    AmiVectorArray const& sx0, AmiVectorArray const& /*sdx0*/
) const {
    int status = CV_SUCCESS;
    sx_ = sx0;
    if (get_sensitivity_method() == SensitivityMethod::forward
        && nplist() > 0) {
        if (get_sens_init_done()) {
            status = CVodeSensReInit(
                solver_memory_.get(),
                static_cast<int>(get_internal_sensitivity_method()),
                sx_.get_nvector_array()
            );
        } else {
            status = CVodeSensInit1(
                solver_memory_.get(), nplist(),
                static_cast<int>(get_internal_sensitivity_method()), fsxdot,
                sx_.get_nvector_array()
            );
            set_sens_init_done();
        }
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensInit1");
}

void CVodeSolver::b_init(
    int const which, realtype const tf, AmiVector const& xB0,
    AmiVector const& /*dxB0*/
) const {
    solver_was_called_B_ = false;
    force_reinit_postprocess_B_ = false;
    xB_ = xB0;
    int status;
    if (get_init_done_b(which)) {
        status
            = CVodeReInitB(solver_memory_.get(), which, tf, xB_.get_nvector());
    } else {
        status = CVodeInitB(
            solver_memory_.get(), which, fxBdot, tf, xB_.get_nvector()
        );
        set_init_done_b(which);
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeInitB");
}

void CVodeSolver::qb_init(int const which, AmiVector const& xQB0) const {
    xQB_ = xQB0;
    int status;
    if (get_quad_init_done_b(which)) {
        status
            = CVodeQuadReInitB(solver_memory_.get(), which, xQB_.get_nvector());
    } else {
        status = CVodeQuadInitB(
            solver_memory_.get(), which, fqBdot, xQB_.get_nvector()
        );
        set_quad_init_done_b(which);
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadInitB");
}

void CVodeSolver::root_init(int ne) const {
    int status = CVodeRootInit(solver_memory_.get(), ne, froot);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeRootInit");
}

void CVodeSolver::set_dense_jac_fn() const {
    int status = CVodeSetJacFn(solver_memory_.get(), fJ);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::set_sparse_jac_fn() const {
    int status = CVodeSetJacFn(solver_memory_.get(), fJSparse);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::set_band_jac_fn() const {
    int status = CVodeSetJacFn(solver_memory_.get(), fJBand);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::set_jac_times_vec_fn() const {
    int status = CVodeSetJacTimes(solver_memory_.get(), nullptr, fJv);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacTimes");
}

void CVodeSolver::set_dense_jac_fn_b(int which) const {
    int status = CVodeSetJacFnB(solver_memory_.get(), which, fJB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::set_sparse_jac_fn_b(int which) const {
    int status = CVodeSetJacFnB(solver_memory_.get(), which, fJSparseB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::set_band_jac_fn_b(int which) const {
    int status = CVodeSetJacFnB(solver_memory_.get(), which, fJBandB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::set_jac_times_vec_fn_b(int which) const {
    int status = CVodeSetJacTimesB(solver_memory_.get(), which, nullptr, fJvB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacTimesB");
}

void CVodeSolver::set_sparse_jac_fn_ss() const {
    int status = CVodeSetJacFn(solver_memory_.get(), fJSparseB_ss);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::apply_max_nonlin_iters() const {
    int status
        = CVodeSetMaxNonlinIters(solver_memory_.get(), get_max_nonlin_iters());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxNonlinIters");
}

void CVodeSolver::apply_max_conv_fails() const {
    int status
        = CVodeSetMaxConvFails(solver_memory_.get(), get_max_conv_fails());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxConvFails");
}

void CVodeSolver::apply_constraints() const {
    Solver::apply_constraints();

    int status = CVodeSetConstraints(
        solver_memory_.get(),
        constraints_.size() > 0 ? constraints_.get_nvector() : nullptr
    );
    if (status != CV_SUCCESS) {
        throw CvodeException(status, "CVodeSetConstraints");
    }
}

void CVodeSolver::apply_max_step_size() const {
    int status = CVodeSetMaxStep(solver_memory_.get(), get_max_step_size());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxStep");
}

Solver* CVodeSolver::clone() const { return new CVodeSolver(*this); }

void CVodeSolver::allocate_solver() const {
    if (!solver_memory_)
        solver_memory_ = std::unique_ptr<void, std::function<void(void*)>>(
            CVodeCreate(static_cast<int>(lmm_), sunctx_),
            [](void* ptr) { CVodeFree(&ptr); }
        );
}

void CVodeSolver::set_ss_tolerances(
    double const rtol, double const atol
) const {
    int status = CVodeSStolerances(solver_memory_.get(), rtol, atol);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSStolerances");
}

void CVodeSolver::set_sens_ss_tolerances(
    double const rtol, double const* atol
) const {
    int status = CVodeSensSStolerances(
        solver_memory_.get(), rtol, const_cast<double*>(atol)
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensSStolerances");
}

void CVodeSolver::set_sens_err_con(bool const error_corr) const {
    int status = CVodeSetSensErrCon(solver_memory_.get(), error_corr);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetSensErrCon");
}

void CVodeSolver::set_quad_err_con_b(int const which, bool const flag) const {
    int status = CVodeSetQuadErrConB(solver_memory_.get(), which, flag);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetQuadErrConB");
}

void CVodeSolver::set_quad_err_con(bool const flag) const {
    int status = CVodeSetQuadErrCon(solver_memory_.get(), flag);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetQuadErrCon");
}

void CVodeSolver::get_root_info(int* rootsfound) const {
    int status = CVodeGetRootInfo(solver_memory_.get(), rootsfound);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetRootInfo");
}

void CVodeSolver::set_linear_solver() const {
    int status = CVodeSetLinearSolver(
        solver_memory_.get(), linear_solver_->get(), linear_solver_->getMatrix()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "setLinearSolver");
}

void CVodeSolver::set_linear_solver_b(int which) const {
    int status = CVodeSetLinearSolverB(
        solver_memory_.get(), which, linear_solver_B_->get(),
        linear_solver_B_->getMatrix()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "setLinearSolverB");
}

void CVodeSolver::set_non_linear_solver() const {
    int status = CVodeSetNonlinearSolver(
        solver_memory_.get(), non_linear_solver_->get()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolver");
}

void CVodeSolver::set_non_linear_solver_sens() const {
    if (get_sensitivity_order() < SensitivityOrder::first)
        return;
    if (get_sensitivity_method() != SensitivityMethod::forward)
        return;

    int status = CV_SUCCESS;

    switch (ism_) {
    case InternalSensitivityMethod::staggered:
        status = CVodeSetNonlinearSolverSensStg(
            solver_memory_.get(), non_linear_solver_sens_->get()
        );
        break;
    case InternalSensitivityMethod::simultaneous:
        status = CVodeSetNonlinearSolverSensSim(
            solver_memory_.get(), non_linear_solver_sens_->get()
        );
        break;
    case InternalSensitivityMethod::staggered1:
        status = CVodeSetNonlinearSolverSensStg1(
            solver_memory_.get(), non_linear_solver_sens_->get()
        );
        break;
    default:
        throw AmiException(
            "Unsupported internal sensitivity method selected: %d", ism_
        );
    }

    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSolver::setNonLinearSolverSens");
}

void CVodeSolver::set_non_linear_solver_b(int const which) const {
    int status = CVodeSetNonlinearSolverB(
        solver_memory_.get(), which, non_linear_solver_B_->get()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolverB");
}

void CVodeSolver::set_user_data() const {
    int status = CVodeSetUserData(solver_memory_.get(), &user_data_);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetUserData");
}

void CVodeSolver::set_user_data_b(int const which) const {
    int status = CVodeSetUserDataB(solver_memory_.get(), which, &user_data_);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetUserDataB");
}

void CVodeSolver::set_max_num_steps(long const mxsteps) const {
    int status = CVodeSetMaxNumSteps(solver_memory_.get(), mxsteps);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxNumSteps");
}

void CVodeSolver::set_stab_lim_det(int const stldet) const {
    int status = CVodeSetStabLimDet(solver_memory_.get(), stldet);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStabLimDet");
}

void CVodeSolver::set_stab_lim_det_b(int const which, int const stldet) const {
    int status = CVodeSetStabLimDetB(solver_memory_.get(), which, stldet);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStabLimDetB");
}

void CVodeSolver::set_id(Model const* /*model*/) const {}

void CVodeSolver::set_suppress_alg(bool const /*flag*/) const {}

void CVodeSolver::reset_state(void* ami_mem, const_N_Vector y0) const {

    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    /* here we force the order in the next step to zero, and update the
     Nordsieck history array, this is largely copied from CVodeReInit with
     explanations from cvodes_impl.h
     */

    /* Set step parameters */

    /* current order */
    cv_mem->cv_q = 1;
    /* L = q + 1 */
    cv_mem->cv_L = 2;
    /* number of steps to wait before updating in q */
    cv_mem->cv_qwait = cv_mem->cv_L;
    /* last successful q value used                */
    cv_mem->cv_qu = 0;
    /* last successful h value used                */
    cv_mem->cv_hu = ZERO;
    /* tolerance scale factor                      */
    cv_mem->cv_tolsf = ONE;

    /* Initialize other integrator optional outputs */

    /* actual initial stepsize                     */
    cv_mem->cv_h0u = ZERO;
    /* step size to be used on the next step       */
    cv_mem->cv_next_h = ZERO;
    /* order to be used on the next step           */
    cv_mem->cv_next_q = 0;

    /* write updated state to Nordsieck history array  */
    N_VScale(ONE, const_cast<N_Vector>(y0), cv_mem->cv_zn[0]);
}

void CVodeSolver::reinit_post_process_f(realtype const tnext) const {
    reInit_post_process(solver_memory_.get(), &t_, &x_, tnext);
    force_reinit_postprocess_F_ = false;
}

void CVodeSolver::reinit_post_process_b(realtype const tnext) const {
    realtype tBret;
    auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
    auto ca_mem = cv_mem->cv_adj_mem;
    auto cvB_mem = ca_mem->cvB_mem;
    // loop over all backward problems
    while (cvB_mem != nullptr) {
        // store current backward problem in ca_mem to make it accessible in
        // adjoint rhs wrapper functions
        ca_mem->ca_bckpbCrt = cvB_mem;
        reInit_post_process(
            static_cast<void*>(cvB_mem->cv_mem), &tBret, &xB_, tnext
        );
        cvB_mem->cv_tout = tBret;
        cvB_mem = cvB_mem->cv_next;
    }
    force_reinit_postprocess_B_ = false;
}

void CVodeSolver::reInit_post_process(
    void* ami_mem, realtype* t, AmiVector* yout, realtype const tout
) const {
    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    auto nst_tmp = cv_mem->cv_nst;
    cv_mem->cv_nst = 0;

    auto status = CVodeSetStopTime(cv_mem, tout);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStopTime");

    status = CVode(ami_mem, tout, yout->get_nvector(), t, CV_ONE_STEP);

    if (status == CV_ROOT_RETURN) {
        auto message
            = std::string("CVode returned a root after reinitialization at t=")
              + std::to_string(*t)
              + ". The initial step-size after the event or "
                "Heaviside function is too small. To fix this, adjust "
                "absolute or relative tolerances!";
        throw CvodeException(status, message.c_str());
    }
    if (status != CV_SUCCESS) {
        std::stringstream msg;
        msg << "tout: " << tout << ", t: " << *t << ".";
        throw CvodeException(status, "reInitPostProcess", msg.str().c_str());
    }

    cv_mem->cv_nst = nst_tmp + 1;
    if (cv_mem->cv_adjMallocDone == SUNTRUE) {
        /* add new step to history array, this is copied from CVodeF */
        auto ca_mem = cv_mem->cv_adj_mem;
        auto dt_mem = ca_mem->dt_mem;

        if (cv_mem->cv_nst % ca_mem->ca_nsteps == 0) {
            /* currently not implemented, we should never get here as we
             limit cv_mem->cv_nst < ca_mem->ca_nsteps, keeping this for
             future regression */
            throw CvodeException(AMICI_ERROR, "reInitPostProcess");
        }

        /* Load next point in dt_mem */
        dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]->t = *t;
        ca_mem->ca_IMstore(cv_mem, dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]);

        /* Set t1 field of the current check point structure
         for the case in which there will be no future
         check points */
        ca_mem->ck_mem->ck_t1 = *t;

        /* tfinal is now set to *tret */
        ca_mem->ca_tfinal = *t;
    }
}

void CVodeSolver::reinit(
    realtype const t0, AmiVector const& yy0, AmiVector const& /*yp0*/
) const {
    auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
    cv_mem->cv_tn = t0;
    if (solver_was_called_F_)
        force_reinit_postprocess_F_ = true;
    x_.copy(yy0);
    reset_state(cv_mem, x_.get_nvector());
}

void CVodeSolver::sens_reinit(
    AmiVectorArray const& yyS0, AmiVectorArray const& /*ypS0*/
) const {
    if (!sens_initialized_) {
        return;
    }

    auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
    /* Initialize znS[0] in the history array */
    for (int is = 0; is < nplist(); is++)
        cv_mem->cv_cvals[is] = ONE;
    if (solver_was_called_F_)
        force_reinit_postprocess_F_ = true;
    sx_.copy(yyS0);
    int status = N_VScaleVectorArray(
        nplist(), cv_mem->cv_cvals, sx_.get_nvector_array(), cv_mem->cv_znS[0]
    );
    if (status != CV_SUCCESS)
        throw CvodeException(CV_VECTOROP_ERR, "CVodeSensReInit");
}

void CVodeSolver::reinit_b(
    int const which, realtype const tB0, AmiVector const& yyB0,
    AmiVector const& /*ypB0*/
) const {
    auto cv_memB = static_cast<CVodeMem>(
        CVodeGetAdjCVodeBmem(solver_memory_.get(), which)
    );
    if (solver_was_called_B_)
        force_reinit_postprocess_B_ = true;
    cv_memB->cv_tn = tB0;
    xB_.copy(yyB0);
    reset_state(cv_memB, xB_.get_nvector());
}

void CVodeSolver::sens_toggle_off() const {
    auto status = CVodeSensToggleOff(solver_memory_.get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensToggleOff");
    /* need to deallocate sensi memory, otherwise can't reenable */
    CVodeSensFree(solver_memory_.get());
    sens_initialized_ = false;
}

void CVodeSolver::reinit_quad_b(int which, AmiVector const& yQB0) const {
    auto cv_memB = static_cast<CVodeMem>(
        CVodeGetAdjCVodeBmem(solver_memory_.get(), which)
    );
    if (solver_was_called_B_)
        force_reinit_postprocess_B_ = true;
    xQB_.copy(yQB0);
    N_VScale(ONE, xQB_.get_nvector(), cv_memB->cv_znQ[0]);
}

void CVodeSolver::set_sens_params(
    realtype const* p, realtype const* pbar, int const* plist
) const {
    int status = CVodeSetSensParams(
        solver_memory_.get(), const_cast<realtype*>(p),
        const_cast<realtype*>(pbar), const_cast<int*>(plist)
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetSensParams");
}

void CVodeSolver::get_dky(realtype t, int k) const {
    int status = CVodeGetDky(solver_memory_.get(), t, k, dky_.get_nvector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetDky");
}

void CVodeSolver::get_sens() const {
    realtype tDummy = 0;
    int status
        = CVodeGetSens(solver_memory_.get(), &tDummy, sx_.get_nvector_array());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetSens");
}

void CVodeSolver::get_sens_dky(realtype const t, int const k) const {
    int status
        = CVodeGetSensDky(solver_memory_.get(), t, k, sx_.get_nvector_array());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetSensDky");
}

void CVodeSolver::get_dky_b(
    realtype const t, int const k, int const which
) const {
    int status = CVodeGetDky(
        CVodeGetAdjCVodeBmem(solver_memory_.get(), which), t, k,
        dky_.get_nvector()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetDkyB");
}

void CVodeSolver::get_quad_b(int which) const {
    realtype tDummy = 0;
    int status = CVodeGetQuadB(
        solver_memory_.get(), which, &tDummy, xQB_.get_nvector()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadB");
}

void CVodeSolver::get_quad(realtype& t) const {
    int status = CVodeGetQuad(solver_memory_.get(), &t, xQ_.get_nvector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuad");
}

void CVodeSolver::get_quad_dky_b(
    realtype const t, int const k, int which
) const {
    int status = CVodeGetQuadDky(
        CVodeGetAdjCVodeBmem(solver_memory_.get(), which), t, k,
        xQB_.get_nvector()
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadDkyB");
}

void CVodeSolver::get_quad_dky(realtype const t, int const k) const {
    int status = CVodeGetQuadDky(solver_memory_.get(), t, k, xQ_.get_nvector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadDky");
}

void CVodeSolver::adj_init() const {
    int status;
    if (get_adj_init_done()) {
        status = CVodeAdjReInit(solver_memory_.get());
    } else {
        status = CVodeAdjInit(
            solver_memory_.get(), static_cast<int>(maxsteps_ + 1),
            static_cast<int>(interp_type_)
        );
        set_adj_init_done();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeAdjInit");
}

void CVodeSolver::quad_init(AmiVector const& xQ0) const {
    int status;
    xQ_.copy(xQ0);
    if (get_quad_init_done()) {
        status = CVodeQuadReInit(
            solver_memory_.get(), const_cast<N_Vector>(xQ0.get_nvector())
        );
    } else {
        status
            = CVodeQuadInit(solver_memory_.get(), fqBdot_ss, xQ_.get_nvector());
        set_quad_init_done();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadInit");
}

void CVodeSolver::allocate_solver_b(int* which) const {
    if (!solver_memory_B_.empty()) {
        *which = 0;
        return;
    }
    int status
        = CVodeCreateB(solver_memory_.get(), static_cast<int>(lmm_), which);
    if (*which + 1 > static_cast<int>(solver_memory_B_.size()))
        solver_memory_B_.resize(*which + 1);
    solver_memory_B_.at(*which)
        = std::unique_ptr<void, std::function<void(void*)>>(
            get_adj_b_mem(solver_memory_.get(), *which), [](void* /*ptr*/) {}
        );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeCreateB");
}

void CVodeSolver::set_ss_tolerances_b(
    int const which, realtype const relTolB, realtype const absTolB
) const {
    int status
        = CVodeSStolerancesB(solver_memory_.get(), which, relTolB, absTolB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSStolerancesB");
}

void CVodeSolver::quad_ss_tolerances_b(
    int const which, realtype const reltolQB, realtype const abstolQB
) const {
    int status = CVodeQuadSStolerancesB(
        solver_memory_.get(), which, reltolQB, abstolQB
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadSStolerancesB");
}

void CVodeSolver::quad_ss_tolerances(
    realtype const reltolQB, realtype const abstolQB
) const {
    int status
        = CVodeQuadSStolerances(solver_memory_.get(), reltolQB, abstolQB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadSStolerances");
}

void CVodeSolver::get_b(int const which) const {
    realtype tDummy = 0;
    int status
        = CVodeGetB(solver_memory_.get(), which, &tDummy, xB_.get_nvector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetB");
}

int CVodeSolver::solve(realtype const tout, int const itask) const {
    if (force_reinit_postprocess_F_)
        reinit_post_process_f(tout);
    int status
        = CVode(solver_memory_.get(), tout, x_.get_nvector(), &t_, itask);
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t_);
    solver_was_called_F_ = true;
    return status;
}

int CVodeSolver::solve_f(
    realtype const tout, int const itask, int* ncheckPtr
) const {
    if (force_reinit_postprocess_F_)
        reinit_post_process_f(tout);
    int status = CVodeF(
        solver_memory_.get(), tout, x_.get_nvector(), &t_, itask, ncheckPtr
    );
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t_);
    solver_was_called_F_ = true;
    return status;
}

void CVodeSolver::solve_b(realtype const tBout, int const itaskB) const {
    if (force_reinit_postprocess_B_)
        reinit_post_process_b(tBout);
    int status = CVodeB(solver_memory_.get(), tBout, itaskB);
    // This does not seem to be documented, but CVodeB may also return
    // CV_TSTOP_RETURN
    // https://github.com/LLNL/sundials/issues/580
    if (status != CV_SUCCESS && status != CV_TSTOP_RETURN) {
        gsl_Expects(status < 0);
        throw IntegrationFailureB(status, tBout);
    }
    solver_was_called_B_ = true;
}

void CVodeSolver::set_max_num_steps_b(
    int const which, long const mxstepsB
) const {
    int status = CVodeSetMaxNumStepsB(solver_memory_.get(), which, mxstepsB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxNumStepsB");
}

void CVodeSolver::diag() const {
    int status = CVDiag(solver_memory_.get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVDiag");
}

void CVodeSolver::diag_b(int const which) const {
    int status = CVDiagB(solver_memory_.get(), which);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVDiagB");
}

void CVodeSolver::get_num_steps(void const* ami_mem, long int* numsteps) const {
    int status = CVodeGetNumSteps(const_cast<void*>(ami_mem), numsteps);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumSteps");
}

void CVodeSolver::get_num_rhs_evals(
    void const* ami_mem, long int* numrhsevals
) const {
    int status = CVodeGetNumRhsEvals(const_cast<void*>(ami_mem), numrhsevals);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumRhsEvals");
}

void CVodeSolver::get_num_err_test_fails(
    void const* ami_mem, long int* numerrtestfails
) const {
    int status
        = CVodeGetNumErrTestFails(const_cast<void*>(ami_mem), numerrtestfails);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumErrTestFails");
}

void CVodeSolver::get_num_non_lin_solv_conv_fails(
    void const* ami_mem, long int* numnonlinsolvconvfails
) const {
    int status = CVodeGetNumNonlinSolvConvFails(
        const_cast<void*>(ami_mem), numnonlinsolvconvfails
    );
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumNonlinSolvConvFails");
}

void CVodeSolver::get_last_order(void const* ami_mem, int* order) const {
    int status = CVodeGetLastOrder(const_cast<void*>(ami_mem), order);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetLastOrder");
}

void* CVodeSolver::get_adj_b_mem(void* ami_mem, int which) const {
    return CVodeGetAdjCVodeBmem(ami_mem, which);
}

void CVodeSolver::calc_ic(realtype const /*tout1*/) const {};

void CVodeSolver::
    calc_ic_b(int const /*which*/, realtype const /*tout1*/) const {};

void CVodeSolver::set_stop_time(realtype const tstop) const {
    int status = CVodeSetStopTime(solver_memory_.get(), tstop);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStopTime");
}

void CVodeSolver::turn_off_root_finding() const {
    int status = CVodeRootInit(solver_memory_.get(), 0, nullptr);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeRootInit");
}

Model const* CVodeSolver::get_model() const {
    if (!solver_memory_)
        throw AmiException(
            "Solver has not been allocated, information is not "
            "available"
        );
    auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());

    auto typed_udata = static_cast<user_data_type*>(cv_mem->cv_user_data);
    Expects(typed_udata);
    return typed_udata->first;
}

/**
 * @brief Jacobian of xdot with respect to states x
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 **/
static int
fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J, void* user_data,
   N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJ(t, x, xdot, J);
    return model->check_finite(J, ModelQuantity::J, t);
}

/**
 * @brief Jacobian of xBdot with respect to adjoint state xB
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 **/
static int
fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector /*tmp1B*/, N_Vector /*tmp2B*/,
    N_Vector /*tmp3B*/) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJB(t, x, xB, xBdot, JB);
    return model->check_finite(gsl::make_span(JB), ModelQuantity::JB, t);
}

/**
 * @brief J in sparse form (for sparse solvers from the SuiteSparse Package)
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 */
static int fJSparse(
    realtype t, N_Vector x, N_Vector /*xdot*/, SUNMatrix J, void* user_data,
    N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/
) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJSparse(t, x, J);
    return model->check_finite(J, ModelQuantity::J, t);
}

/**
 * @brief JB in sparse form (for sparse solvers from the SuiteSparse Package)
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 */
static int fJSparseB(
    realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector /*tmp1B*/, N_Vector /*tmp2B*/, N_Vector /*tmp3B*/
) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJSparseB(t, x, xB, xBdot, JB);
    return model->check_finite(gsl::make_span(JB), ModelQuantity::JB, t);
}

/**
 * @brief J in banded form (for banded solvers)
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 */
static int fJBand(
    realtype t, N_Vector x, N_Vector xdot, SUNMatrix J, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
) {
    return fJ(t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}

/**
 * @brief JB in banded form (for banded solvers)
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 */
static int fJBandB(
    realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB,
    void* user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B
) {
    return fJB(t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

/**
 * @brief Matrix vector product of J with a vector v (for iterative solvers)
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param v Vector with which the Jacobian is multiplied
 * @param Jv Vector to which the Jacobian vector product will be
 *written
 * @param user_data object with user input
 * @param tmp temporary storage vector
 * @return status flag indicating successful execution
 **/
static int
fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector /*xdot*/,
    void* user_data, N_Vector /*tmp*/) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJv(v, Jv, t, x);
    return model->check_finite(gsl::make_span(Jv), ModelQuantity::Jv, t);
}

/**
 * @brief Matrix vector product of JB with a vector v (for iterative solvers)
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param vB Vector with which the Jacobian is multiplied
 * @param JvB Vector to which the Jacobian vector product will be
 *written
 * @param user_data object with user input
 * @param tmpB temporary storage vector
 * @return status flag indicating successful execution
 **/
static int fJvB(
    N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB,
    N_Vector /*xBdot*/, void* user_data, N_Vector /*tmpB*/
) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJvB(vB, JvB, t, x, xB);
    return model->check_finite(gsl::make_span(JvB), ModelQuantity::JvB, t);
}

/**
 * @brief Event trigger function for events
 * @param t timepoint
 * @param x Vector with the states
 * @param root array with root function values
 * @param user_data object with user input
 * @return status flag indicating successful execution
 */
static int froot(realtype t, N_Vector x, realtype* root, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    if (model->ne != model->ne_solver) {
        // temporary buffer to store all root function values, not only the ones
        // tracked by the solver
        thread_local static std::vector<realtype> root_buffer(model->ne, 0.0);
        root_buffer.resize(model->ne);
        model->froot(t, x, root_buffer);
        std::copy_n(root_buffer.begin(), model->ne_solver, root);
    } else {
        model->froot(t, x, gsl::make_span<realtype>(root, model->ne_solver));
    }
    return model->check_finite(
        gsl::make_span<realtype>(root, model->ne_solver), ModelQuantity::root, t
    );
}

/**
 * @brief residual function of the ODE
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param user_data object with user input
 * @return status flag indicating successful execution
 */
static int fxdot(realtype t, N_Vector x, N_Vector xdot, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);
    auto solver = dynamic_cast<CVodeSolver const*>(typed_udata->second);
    Expects(model);

    if (solver->time_exceeded(500)) {
        return AMICI_MAX_TIME_EXCEEDED;
    }


    model->fxdot(t, x, xdot);
    return model->check_finite(gsl::make_span(xdot), ModelQuantity::xdot, t);
}

/**
 * @brief Right hand side of differential equation for adjoint state xB
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param user_data object with user input
 * @return status flag indicating successful execution
 */
static int
fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);
    auto solver = dynamic_cast<CVodeSolver const*>(typed_udata->second);
    Expects(model);

    if (solver->time_exceeded(500)) {
        return AMICI_MAX_TIME_EXCEEDED;
    }

    model->fxBdot(t, x, xB, xBdot);
    return model->check_finite(gsl::make_span(xBdot), ModelQuantity::xBdot, t);
}

/**
 * @brief Right hand side of integral equation for quadrature states qB
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param qBdot Vector with the adjoint quadrature right hand side
 * @param user_data pointer to temp data object
 * @return status flag indicating successful execution
 */
static int
fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fqBdot(t, x, xB, qBdot);
    return model->check_finite(gsl::make_span(qBdot), ModelQuantity::qBdot, t);
}

/**
 * @brief Right hand side of differential equation for adjoint state xB
 * when simulating in steadystate mode
 * @param t timepoint
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param user_data object with user input
 * @return status flag indicating successful execution
 */
static int fxBdot_ss(realtype t, N_Vector xB, N_Vector xBdot, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fxBdot_ss(t, xB, xBdot);
    return model->check_finite(
        gsl::make_span(xBdot), ModelQuantity::xBdot_ss, t
    );
}

/**
 * @brief Right hand side of integral equation for quadrature states qB
 * when simulating in steadystate mode
 * @param t timepoint
 * @param xB Vector with the adjoint states
 * @param qBdot Vector with the adjoint quadrature right hand side
 * @param user_data pointer to temp data object
 * @return status flag indicating successful execution
 */
static int fqBdot_ss(realtype t, N_Vector xB, N_Vector qBdot, void* user_data) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fqBdot_ss(t, xB, qBdot);
    return model->check_finite(
        gsl::make_span(qBdot), ModelQuantity::qBdot_ss, t
    );
}

/**
 * @brief JB in sparse form for steady state case
 * @param t timepoint
 * @param x Vector with the states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 */
static int fJSparseB_ss(
    realtype t, N_Vector /*x*/, N_Vector xBdot, SUNMatrix JB, void* user_data,
    N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/
) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fJSparseB_ss(JB);
    return model->check_finite(
        gsl::make_span(xBdot), ModelQuantity::JSparseB_ss, t
    );
}

/**
 * @brief Right hand side of differential equation for state sensitivities sx
 * @param Ns number of parameters
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param ip parameter index
 * @param sx Vector with the state sensitivities
 * @param sxdot Vector with the sensitivity right hand side
 * @param user_data object with user input
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 */
static int fsxdot(
    int /*Ns*/, realtype t, N_Vector x, N_Vector /*xdot*/, int ip, N_Vector sx,
    N_Vector sxdot, void* user_data, N_Vector /*tmp1*/, N_Vector /*tmp2*/
) {
    auto typed_udata = static_cast<CVodeSolver::user_data_type*>(user_data);
    Expects(typed_udata);
    auto model = dynamic_cast<Model_ODE*>(typed_udata->first);
    Expects(model);

    model->fsxdot(t, x, ip, sx, sxdot);
    return model->check_finite(gsl::make_span(sxdot), ModelQuantity::sxdot, t);
}

bool operator==(CVodeSolver const& a, CVodeSolver const& b) {
    return static_cast<Solver const&>(a) == static_cast<Solver const&>(b);
}

} // namespace amici
