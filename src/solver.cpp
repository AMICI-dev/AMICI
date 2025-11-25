#include "amici/solver.h"

#include "amici/amici.h"
#include "amici/exception.h"
#include "amici/model.h"

#include <sundials/sundials_context.h>
#include <gsl/gsl-lite.hpp>

#include <cstdio>
#include <filesystem>
#include <memory>

namespace amici {
using std::isnan;

/* Error handler passed to SUNDIALS. */
void wrapErrHandlerFn(
    [[maybe_unused]] int const line, char const* func, char const* file,
    char const* msg, SUNErrCode const err_code, void* err_user_data,
    [[maybe_unused]] SUNContext sunctx
) {
    constexpr int BUF_SIZE = 250;

    char msg_buffer[BUF_SIZE];
    char id_buffer[BUF_SIZE];
    static_assert(std::is_same_v<SUNErrCode, int>, "Must update format string");
    // for debug builds, include full file path and line numbers
#ifndef NDEBUG
    snprintf(msg_buffer, BUF_SIZE, "%s:%d: %s (%d)", file, line, msg, err_code);
#else
    snprintf(msg_buffer, BUF_SIZE, "%s", msg);
#endif
    // we need a matlab-compatible message ID
    // i.e. colon separated and only  [A-Za-z0-9_]
    // see https://mathworks.com/help/matlab/ref/mexception.html
    std::filesystem::path const path(file);
    auto file_stem = path.stem().string();

    // Error code to string. Remove 'AMICI_' prefix.
    auto err_code_str = simulation_status_to_str(err_code);
    constexpr std::string_view err_code_prefix = "AMICI_";
    if (err_code_str.substr(0, err_code_prefix.size()) == err_code_prefix)
        err_code_str = err_code_str.substr(err_code_prefix.size());

    snprintf(
        id_buffer, BUF_SIZE, "%s:%s:%s", file_stem.c_str(), func,
        err_code_str.c_str()
    );

    if (!err_user_data) {
        throw std::runtime_error("eh_data unset");
    }

    auto solver = static_cast<Solver const*>(err_user_data);
    if (solver->get_logger())
        solver->get_logger()->log(LogSeverity::debug, id_buffer, msg_buffer);
}

Solver::Solver(Solver const& other)
    : ism_(other.ism_)
    , lmm_(other.lmm_)
    , iter_(other.iter_)
    , interp_type_(other.interp_type_)
    , maxsteps_(other.maxsteps_)
    , maxtime_(other.maxtime_)
    , simulation_timer_(other.simulation_timer_)
    , constraints_(other.constraints_)
    , sensi_meth_(other.sensi_meth_)
    , sensi_meth_preeq_(other.sensi_meth_preeq_)
    , stldet_(other.stldet_)
    , ordering_(other.ordering_)
    , newton_maxsteps_(other.newton_maxsteps_)
    , newton_damping_factor_mode_(other.newton_damping_factor_mode_)
    , newton_damping_factor_lower_bound_(
          other.newton_damping_factor_lower_bound_
      )
    , linsol_(other.linsol_)
    , atol_(other.atol_)
    , rtol_(other.rtol_)
    , atol_fsa_(other.atol_fsa_)
    , rtol_fsa_(other.rtol_fsa_)
    , atolB_(other.atolB_)
    , rtolB_(other.rtolB_)
    , quad_atol_(other.quad_atol_)
    , quad_rtol_(other.quad_rtol_)
    , ss_tol_factor_(other.ss_tol_factor_)
    , ss_atol_(other.ss_atol_)
    , ss_rtol_(other.ss_rtol_)
    , ss_tol_sensi_factor_(other.ss_tol_sensi_factor_)
    , ss_atol_sensi_(other.ss_atol_sensi_)
    , ss_rtol_sensi_(other.ss_rtol_sensi_)
    , rdata_mode_(other.rdata_mode_)
    , newton_step_steadystate_conv_(other.newton_step_steadystate_conv_)
    , check_sensi_steadystate_conv_(other.check_sensi_steadystate_conv_)
    , max_nonlin_iters_(other.max_nonlin_iters_)
    , max_conv_fails_(other.max_conv_fails_)
    , max_step_size_(other.max_step_size_)
    , maxstepsB_(other.maxstepsB_)
    , sensi_(other.sensi_) {
    // update to our own context
    constraints_.set_ctx(sunctx_);
}

SUNContext Solver::get_sun_context() const { return sunctx_; }

void Solver::apply_max_num_steps() const {
    // set remaining steps, setMaxNumSteps only applies to a single call of
    // solve
    long int cursteps;
    get_num_steps(solver_memory_.get(), &cursteps);
    if (maxsteps_ <= cursteps)
        throw AmiException(
            "Reached maximum number of steps %ld before reaching "
            "tout at t=%g.",
            maxsteps_, t_
        );
    set_max_num_steps(maxsteps_ - cursteps);
}

void Solver::apply_max_num_steps_B() const {
    // set remaining steps, setMaxNumSteps only applies to a single call of
    // solve
    long int curstepsB;
    auto maxstepsB = (maxstepsB_ == 0) ? maxsteps_ * 100 : maxstepsB_;
    for (int i_mem_b = 0; i_mem_b < (int)solver_memory_B_.size(); ++i_mem_b) {
        if (solver_memory_B_.at(i_mem_b)) {
            get_num_steps(solver_memory_B_.at(i_mem_b).get(), &curstepsB);
            if (maxstepsB <= curstepsB)
                throw AmiException(
                    "Reached maximum number of steps %ld before "
                    "reaching tout at t=%g in backward "
                    "problem %i.",
                    maxstepsB_, t_, i_mem_b
                );
            set_max_num_steps_b(i_mem_b, maxstepsB - curstepsB);
        }
    }
}

int Solver::run(realtype const tout) const {
    gsl_ExpectsDebug(std::isfinite(tout));

    set_stop_time(tout);
    CpuTimer const cpu_timer;
    int status = AMICI_SUCCESS;

    apply_max_num_steps();
    if (nx() > 0) {
        if (get_adj_init_done()) {
            status = solve_f(tout, AMICI_NORMAL, &ncheckPtr_);
        } else {
            status = solve(tout, AMICI_NORMAL);
        }
    } else {
        t_ = tout;
    }
    cpu_time_ += cpu_timer.elapsed_milliseconds();
    return status;
}

int Solver::step(realtype const tout) const {
    gsl_ExpectsDebug(std::isfinite(tout));

    int status = AMICI_SUCCESS;

    apply_max_num_steps();
    if (nx() > 0) {
        if (get_adj_init_done()) {
            status = solve_f(tout, AMICI_ONE_STEP, &ncheckPtr_);
        } else {
            status = solve(tout, AMICI_ONE_STEP);
        }
    } else {
        t_ = tout;
    }
    return status;
}

void Solver::run_b(realtype const tout) const {
    gsl_ExpectsDebug(std::isfinite(tout));

    CpuTimer cpu_timer;

    apply_max_num_steps_B();
    if (nx() > 0) {
        solve_b(tout, AMICI_NORMAL);
    }
    cpu_time_b_ += cpu_timer.elapsed_milliseconds();
    t_ = tout;
}

void Solver::setup(
    realtype const t0, Model* model, AmiVector const& x0, AmiVector const& dx0,
    AmiVectorArray const& sx0, AmiVectorArray const& sdx0
) const {
    if (nx() != model->nx_solver || nplist() != model->nplist()
        || nquad() != model->nJ * model->nplist()) {
        reset_mutable_memory(
            model->nx_solver, model->nplist(), model->nJ * model->nplist()
        );
    }
    /* Create solver memory object if necessary */
    allocate_solver();
    if (!solver_memory_)
        throw AmiException("Failed to allocated solver memory!");

    /* Initialize CVodes/IDAs solver*/
    init(t0, x0, dx0);

    /* Clear diagnosis storage */
    reset_diagnosis();

    /* Apply stored tolerances to sundials solver */
    apply_tolerances();

    /* Set optional inputs */
    set_err_handler_fn();
    /* Attaches userdata */
    user_data_ = std::make_pair(model, this);
    set_user_data();
    /* activates stability limit detection */
    set_stab_lim_det(stldet_);

    root_init(model->ne_solver);

    if (nx() == 0)
        return;

    initialize_linear_solver(model);
    initialize_non_linear_solver();

    if (sensi_ >= SensitivityOrder::first
        && sensi_meth_ > SensitivityMethod::none && model->nx_solver > 0) {
        auto const plist = model->get_parameter_list();
        sens_init_1(sx0, sdx0);
        if (sensi_meth_ == SensitivityMethod::forward && !plist.empty()) {
            /* Set sensitivity analysis optional inputs */
            auto const par = model->get_unscaled_parameters();

            /* Activate sensitivity calculations  and apply tolerances */
            initialize_non_linear_solver_sens(model);
            set_sens_params(par.data(), nullptr, plist.data());
            apply_tolerances_fsa();
        } else {
            /* Allocate space for the adjoint computation */
            adj_init();
        }
    }

    set_id(model);
    set_suppress_alg(true);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    if (model->nt() > 1)
        calc_ic(model->get_timepoint(1));

    apply_max_nonlin_iters();
    apply_max_conv_fails();
    apply_max_step_size();

    cpu_time_ = 0.0;
    cpu_time_b_ = 0.0;

    apply_constraints();
}

void Solver::setup_b(
    int* which, realtype const tf, Model* model, AmiVector const& xB0,
    AmiVector const& dxB0, AmiVector const& xQB0
) const {
    if (!solver_memory_)
        throw AmiException(
            "Solver for the forward problem must be setup first"
        );

    /* allocate memory for the backward problem */
    allocate_solver_b(which);

    /* initialise states */
    b_init(*which, tf, xB0, dxB0);

    /* Attach user data */
    set_user_data_b(*which);

    if (nx() == 0)
        return;

    initialize_linear_solver_b(model, *which);
    initialize_non_linear_solver_b(*which);

    /* Initialise quadrature calculation */
    qb_init(*which, xQB0);

    apply_tolerances_asa(*which);
    apply_quad_tolerances_asa(*which);

    set_stab_lim_det_b(*which, stldet_);
}

void Solver::setup_steady_state(
    realtype const t0, Model* model, AmiVector const& x0, AmiVector const& dx0,
    AmiVector const& xB0, AmiVector const& dxB0, AmiVector const& xQ0
) const {
    /* Initialize CVodes/IDAs solver with steadystate RHS function */
    init_steady_state(t0, x0, dx0);

    /* Allocate space for forward quadratures */
    quad_init(xQ0);

    /* Apply tolerances */
    apply_quad_tolerances();

    /* Check linear solver (works only with KLU atm) */
    if (linsol_ != LinearSolver::KLU)
        throw AmiException(
            "Backward steady state computation via integration "
            "is currently only implemented for KLU linear solver"
        );
    /* Set Jacobian function and initialize values */
    set_sparse_jac_fn_ss();
    model->write_steady_state_JB(t0, 0, x0, dx0, xB0, dxB0, xB0);
}

void Solver::update_and_reinit_states_and_sensitivities(Model* model) const {
    model->reinitialize(
        t_, x_, sx_, get_sensitivity_order() >= SensitivityOrder::first
    );

    reinit(t_, x_, dx_);
    if (get_sensitivity_order() >= SensitivityOrder::first
        && get_sensitivity_method() == SensitivityMethod::forward) {
        sens_reinit(sx_, sdx_);
    }
}

void Solver::reset_diagnosis() const {
    ns_.clear();
    nrhs_.clear();
    netf_.clear();
    nnlscf_.clear();
    order_.clear();

    nsB_.clear();
    nrhsB_.clear();
    netfB_.clear();
    nnlscfB_.clear();
}

void Solver::store_diagnosis() const {
    if (!solver_was_called_F_ || !solver_memory_) {
        ns_.push_back(0);
        nrhs_.push_back(0);
        netf_.push_back(0);
        nnlscf_.push_back(0);
        order_.push_back(0);
        return;
    }

    long int lnumber;
    get_num_steps(solver_memory_.get(), &lnumber);
    ns_.push_back(gsl::narrow<int>(lnumber));

    get_num_rhs_evals(solver_memory_.get(), &lnumber);
    nrhs_.push_back(gsl::narrow<int>(lnumber));

    get_num_err_test_fails(solver_memory_.get(), &lnumber);
    netf_.push_back(gsl::narrow<int>(lnumber));

    get_num_non_lin_solv_conv_fails(solver_memory_.get(), &lnumber);
    nnlscf_.push_back(gsl::narrow<int>(lnumber));

    int number;
    get_last_order(solver_memory_.get(), &number);
    order_.push_back(number);
}

void Solver::store_diagnosis_b(int const which) const {
    if (!solver_was_called_B_ || !solver_memory_B_.at(which)) {
        nsB_.push_back(0);
        nrhsB_.push_back(0);
        netfB_.push_back(0);
        nnlscfB_.push_back(0);
        return;
    }

    long int number;
    get_num_steps(solver_memory_B_.at(which).get(), &number);
    nsB_.push_back(gsl::narrow<int>(number));

    get_num_rhs_evals(solver_memory_B_.at(which).get(), &number);
    nrhsB_.push_back(gsl::narrow<int>(number));

    get_num_err_test_fails(solver_memory_B_.at(which).get(), &number);
    netfB_.push_back(gsl::narrow<int>(number));

    get_num_non_lin_solv_conv_fails(solver_memory_B_.at(which).get(), &number);
    nnlscfB_.push_back(gsl::narrow<int>(number));
}

void Solver::initialize_linear_solver(Model const* model) const {
    switch (linsol_) {

        /* DIRECT SOLVERS */

    case LinearSolver::dense:
        linear_solver_ = std::make_unique<SUNLinSolDense>(x_);
        set_linear_solver();
        set_dense_jac_fn();
        break;

    case LinearSolver::band:
        linear_solver_
            = std::make_unique<SUNLinSolBand>(x_, model->ubw, model->lbw);
        set_linear_solver();
        set_band_jac_fn();
        break;

    case LinearSolver::LAPACKDense:
    case LinearSolver::LAPACKBand:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::diag:
        diag();
        set_dense_jac_fn();
        break;

        /* ITERATIVE SOLVERS */

    case LinearSolver::SPGMR:
        linear_solver_ = std::make_unique<SUNLinSolSPGMR>(x_);
        set_linear_solver();
        set_jac_times_vec_fn();
        break;

    case LinearSolver::SPBCG:
        linear_solver_ = std::make_unique<SUNLinSolSPBCGS>(x_);
        set_linear_solver();
        set_jac_times_vec_fn();
        break;

    case LinearSolver::SPTFQMR:
        linear_solver_ = std::make_unique<SUNLinSolSPTFQMR>(x_);
        set_linear_solver();
        set_jac_times_vec_fn();
        break;

        /* SPARSE SOLVERS */

    case LinearSolver::KLU:
        linear_solver_ = std::make_unique<SUNLinSolKLU>(
            x_, model->nnz, CSC_MAT,
            static_cast<SUNLinSolKLU::StateOrdering>(get_state_ordering())
        );
        set_linear_solver();
        set_sparse_jac_fn();
        break;

#ifdef SUNDIALS_SUPERLUMT
    case LinearSolver::SuperLUMT:
        // TODO state ordering
        linearSolver = std::make_unique<SUNLinSolSuperLUMT>(
            *x, model->nnz, CSC_MAT,
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(getStateOrdering())
        );

        setLinearSolver();
        setSparseJacFn();
        break;
#endif
    default:
        throw AmiException(
            "Invalid choice of solver: %d", static_cast<int>(linsol_)
        );
    }
}

void Solver::initialize_non_linear_solver() const {
    switch (iter_) {
    case NonlinearSolverIteration::newton:
        non_linear_solver_
            = std::make_unique<SUNNonLinSolNewton>(x_.get_nvector());
        break;
    case NonlinearSolverIteration::fixedpoint:
        non_linear_solver_
            = std::make_unique<SUNNonLinSolFixedPoint>(x_.get_nvector());
        break;
    default:
        throw AmiException(
            "Invalid non-linear solver specified (%d).", static_cast<int>(iter_)
        );
    }

    set_non_linear_solver();
}

void Solver::initialize_linear_solver_b(
    Model const* model, int const which
) const {
    switch (linsol_) {
    /* DIRECT SOLVERS */
    case LinearSolver::dense:
        linear_solver_B_ = std::make_unique<SUNLinSolDense>(xB_);
        set_linear_solver_b(which);
        set_dense_jac_fn_b(which);
        break;

    case LinearSolver::band:
        linear_solver_B_
            = std::make_unique<SUNLinSolBand>(xB_, model->ubw, model->lbw);
        set_linear_solver_b(which);
        set_band_jac_fn_b(which);
        break;

    case LinearSolver::LAPACKDense:
    case LinearSolver::LAPACKBand:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::diag:
        diag_b(which);
        set_dense_jac_fn_b(which);
        break;

        /* ITERATIVE SOLVERS */

    case LinearSolver::SPGMR:
        linear_solver_B_ = std::make_unique<SUNLinSolSPGMR>(xB_);
        set_linear_solver_b(which);
        set_jac_times_vec_fn_b(which);
        break;

    case LinearSolver::SPBCG:
        linear_solver_B_ = std::make_unique<SUNLinSolSPBCGS>(xB_);
        set_linear_solver_b(which);
        set_jac_times_vec_fn_b(which);
        break;

    case LinearSolver::SPTFQMR:
        linear_solver_B_ = std::make_unique<SUNLinSolSPTFQMR>(xB_);
        set_linear_solver_b(which);
        set_jac_times_vec_fn_b(which);
        break;

        /* SPARSE SOLVERS */

    case LinearSolver::KLU:
        linear_solver_B_ = std::make_unique<SUNLinSolKLU>(
            xB_, model->nnz, CSC_MAT,
            static_cast<SUNLinSolKLU::StateOrdering>(get_state_ordering())
        );
        set_linear_solver_b(which);
        set_sparse_jac_fn_b(which);
        break;
#ifdef SUNDIALS_SUPERLUMT
    case LinearSolver::SuperLUMT:
        linearSolverB = std::make_unique<SUNLinSolSuperLUMT>(
            *xB, model->nnz, CSC_MAT,
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(getStateOrdering())
        );
        setLinearSolverB(which);
        setSparseJacFnB(which);
        break;
#endif
    default:
        throw AmiException(
            "Invalid choice of solver: %d", static_cast<int>(linsol_)
        );
    }
}

void Solver::initialize_non_linear_solver_b(int const which) const {
    switch (iter_) {
    case NonlinearSolverIteration::newton:
        non_linear_solver_B_
            = std::make_unique<SUNNonLinSolNewton>(xB_.get_nvector());
        break;
    case NonlinearSolverIteration::fixedpoint:
        non_linear_solver_B_
            = std::make_unique<SUNNonLinSolFixedPoint>(xB_.get_nvector());
        break;
    default:
        throw AmiException(
            "Invalid non-linear solver specified (%d).", static_cast<int>(iter_)
        );
    }

    set_non_linear_solver_b(which);
}

bool operator==(Solver const& a, Solver const& b) {
    if (typeid(a) != typeid(b))
        return false;

    return (a.interp_type_ == b.interp_type_) && (a.lmm_ == b.lmm_)
           && (a.iter_ == b.iter_) && (a.stldet_ == b.stldet_)
           && (a.ordering_ == b.ordering_)
           && (a.newton_maxsteps_ == b.newton_maxsteps_)
           && (a.newton_damping_factor_mode_ == b.newton_damping_factor_mode_)
           && (a.newton_damping_factor_lower_bound_
               == b.newton_damping_factor_lower_bound_)
           && (a.ism_ == b.ism_) && (a.linsol_ == b.linsol_)
           && (a.atol_ == b.atol_) && (a.rtol_ == b.rtol_)
           && (a.maxsteps_ == b.maxsteps_) && (a.maxstepsB_ == b.maxstepsB_)
           && (a.quad_atol_ == b.quad_atol_) && (a.quad_rtol_ == b.quad_rtol_)
           && (a.maxtime_ == b.maxtime_)
           && (a.get_absolute_tolerance_steady_state()
               == b.get_absolute_tolerance_steady_state())
           && (a.get_relative_tolerance_steady_state()
               == b.get_relative_tolerance_steady_state())
           && (a.get_absolute_tolerance_steady_state_sensi()
               == b.get_absolute_tolerance_steady_state_sensi())
           && (a.get_relative_tolerance_steady_state_sensi()
               == b.get_relative_tolerance_steady_state_sensi())
           && (a.rtol_fsa_ == b.rtol_fsa_
               || (isnan(a.rtol_fsa_) && isnan(b.rtol_fsa_)))
           && (a.atol_fsa_ == b.atol_fsa_
               || (isnan(a.atol_fsa_) && isnan(b.atol_fsa_)))
           && (a.rtolB_ == b.rtolB_ || (isnan(a.rtolB_) && isnan(b.rtolB_)))
           && (a.atolB_ == b.atolB_ || (isnan(a.atolB_) && isnan(b.atolB_)))
           && (a.sensi_ == b.sensi_) && (a.sensi_meth_ == b.sensi_meth_)
           && (a.newton_step_steadystate_conv_
               == b.newton_step_steadystate_conv_)
           && (a.check_sensi_steadystate_conv_
               == b.check_sensi_steadystate_conv_)
           && (a.rdata_mode_ == b.rdata_mode_)
           && (a.max_conv_fails_ == b.max_conv_fails_)
           && (a.max_nonlin_iters_ == b.max_nonlin_iters_)
           && (a.max_step_size_ == b.max_step_size_)
           && (a.constraints_.get_vector() == b.constraints_.get_vector());
}

void Solver::apply_tolerances() const {
    if (!get_init_done())
        throw AmiException(
            "Solver instance was not yet set up, the "
            "tolerances cannot be applied yet!"
        );

    set_ss_tolerances(rtol_, atol_);
}

void Solver::apply_tolerances_fsa() const {
    if (!get_init_done())
        throw AmiException(
            "Solver instance was not yet set up, the "
            "tolerances cannot be applied yet!"
        );

    if (sensi_ < SensitivityOrder::first)
        return;

    if (nplist()) {
        std::vector<realtype> atols(nplist(), get_absolute_tolerance_fsa());
        set_sens_ss_tolerances(get_relative_tolerance_fsa(), atols.data());
        set_sens_err_con(true);
    }
}

void Solver::apply_tolerances_asa(int const which) const {
    if (!get_adj_init_done())
        throw AmiException(
            "Adjoint solver instance was not yet set up, the "
            "tolerances cannot be applied yet!"
        );

    if (sensi_ < SensitivityOrder::first)
        return;

    /* specify integration tolerances for backward problem */
    set_ss_tolerances_b(
        which, get_relative_tolerance_b(), get_absolute_tolerance_b()
    );
}

void Solver::apply_quad_tolerances_asa(int const which) const {
    if (!get_adj_init_done())
        throw AmiException(
            "Adjoint solver instance was not yet set up, the "
            "tolerances cannot be applied yet!"
        );

    if (sensi_ < SensitivityOrder::first)
        return;

    realtype const quad_rtol = isnan(quad_rtol_) ? rtol_ : quad_rtol_;
    realtype const quad_atol = isnan(quad_atol_) ? atol_ : quad_atol_;

    /* Enable Quadrature Error Control */
    set_quad_err_con_b(which, !std::isinf(quad_atol) && !std::isinf(quad_rtol));

    quad_ss_tolerances_b(which, quad_rtol, quad_atol);
}

void Solver::apply_quad_tolerances() const {
    if (!get_quad_init_done())
        throw AmiException(
            "Quadratures were not initialized, the "
            "tolerances cannot be applied yet!"
        );

    if (sensi_ < SensitivityOrder::first)
        return;

    realtype quad_rtolF = isnan(quad_rtol_) ? rtol_ : quad_rtol_;
    realtype quad_atolF = isnan(quad_atol_) ? atol_ : quad_atol_;

    /* Enable Quadrature Error Control */
    set_quad_err_con(!std::isinf(quad_atolF) && !std::isinf(quad_rtolF));

    quad_ss_tolerances(quad_rtolF, quad_atolF);
}

void Solver::apply_sensitivity_tolerances() const {
    if (sensi_ < SensitivityOrder::first)
        return;

    if (sensi_meth_ == SensitivityMethod::forward)
        apply_tolerances_fsa();
    else if (sensi_meth_ == SensitivityMethod::adjoint && get_adj_init_done()) {
        for (int iMem = 0; iMem < (int)solver_memory_B_.size(); ++iMem)
            apply_tolerances_asa(iMem);
    }
}

void Solver::apply_constraints() const {
    if (constraints_.size() != 0
        && gsl::narrow<int>(constraints_.size()) != nx()) {
        throw std::invalid_argument(
            "Constraints must have the same size as the state vector."
        );
    }
}

SensitivityMethod Solver::get_sensitivity_method() const { return sensi_meth_; }

SensitivityMethod Solver::get_sensitivity_method_pre_equilibration() const {
    return sensi_meth_preeq_;
}

void Solver::set_sensitivity_method(SensitivityMethod const sensi_meth) {
    check_sensitivity_method(sensi_meth, false);
    this->sensi_meth_ = sensi_meth;
}

void Solver::set_sensitivity_method_pre_equilibration(
    SensitivityMethod const sensi_meth_preeq
) {
    check_sensitivity_method(sensi_meth_preeq, true);
    sensi_meth_preeq_ = sensi_meth_preeq;
}

void Solver::check_sensitivity_method(
    SensitivityMethod const sensi_meth, bool const preequilibration
) const {
    if (rdata_mode_ == RDataReporting::residuals
        && sensi_meth == SensitivityMethod::adjoint)
        throw AmiException(
            "Adjoint Sensitivity Analysis is not compatible with"
            " only reporting residuals!"
        );
    if (!preequilibration && sensi_meth != sensi_meth_)
        reset_mutable_memory(nx(), nplist(), nquad());
}

void Solver::set_max_nonlin_iters(int const max_nonlin_iters) {
    if (max_nonlin_iters < 0)
        throw AmiException("max_nonlin_iters must be a non-negative number");

    max_nonlin_iters_ = max_nonlin_iters;
}

int Solver::get_max_nonlin_iters() const { return max_nonlin_iters_; }

void Solver::set_max_conv_fails(int const max_conv_fails) {
    if (max_conv_fails < 0)
        throw AmiException("max_conv_fails must be a non-negative number");

    max_conv_fails_ = max_conv_fails;
}

int Solver::get_max_conv_fails() const { return max_conv_fails_; }

void Solver::set_constraints(std::vector<realtype> const& constraints) {
    auto any_constraint
        = std::ranges::any_of(constraints, [](bool x) { return x != 0.0; });

    if (!any_constraint) {
        // all-0 must be converted to empty, otherwise sundials will fail
        constraints_ = AmiVector(0, sunctx_);
        return;
    }

    constraints_ = AmiVector(constraints, sunctx_);
}

void Solver::set_max_step_size(realtype max_step_size) {
    if (max_step_size < 0)
        throw AmiException("max_step_size must be non-negative.");
    max_step_size_ = max_step_size;
}

realtype Solver::get_max_step_size() const { return max_step_size_; }

int Solver::get_newton_max_steps() const { return newton_maxsteps_; }

void Solver::set_newton_max_steps(int const newton_maxsteps) {
    if (newton_maxsteps < 0)
        throw AmiException("newton_maxsteps must be a non-negative number");
    newton_maxsteps_ = newton_maxsteps;
}

NewtonDampingFactorMode Solver::get_newton_damping_factor_mode() const {
    return newton_damping_factor_mode_;
}

void Solver::set_newton_damping_factor_mode(
    NewtonDampingFactorMode dampingFactorMode
) {
    newton_damping_factor_mode_ = dampingFactorMode;
}

double Solver::get_newton_damping_factor_lower_bound() const {
    return newton_damping_factor_lower_bound_;
}

void Solver::set_newton_damping_factor_lower_bound(
    double const dampingFactorLowerBound
) {
    newton_damping_factor_lower_bound_ = dampingFactorLowerBound;
}

SensitivityOrder Solver::get_sensitivity_order() const { return sensi_; }

void Solver::set_sensitivity_order(SensitivityOrder const sensi) {
    if (sensi_ != sensi)
        reset_mutable_memory(nx(), nplist(), nquad());
    sensi_ = sensi;

    if (get_init_done())
        apply_sensitivity_tolerances();
}

double Solver::get_relative_tolerance() const {
    return static_cast<double>(rtol_);
}

void Solver::set_relative_tolerance(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtol_ = static_cast<realtype>(rtol);

    if (get_init_done()) {
        apply_tolerances();
        apply_sensitivity_tolerances();
    }
}

double Solver::get_absolute_tolerance() const {
    return static_cast<double>(atol_);
}

void Solver::set_absolute_tolerance(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    atol_ = static_cast<realtype>(atol);

    if (get_init_done()) {
        apply_tolerances();
        apply_sensitivity_tolerances();
    }
}

double Solver::get_relative_tolerance_fsa() const {
    return static_cast<double>(isnan(rtol_fsa_) ? rtol_ : rtol_fsa_);
}

void Solver::set_relative_tolerance_fsa(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtol_fsa_ = static_cast<realtype>(rtol);

    if (get_init_done()) {
        apply_sensitivity_tolerances();
    }
}

double Solver::get_absolute_tolerance_fsa() const {
    return static_cast<double>(isnan(atol_fsa_) ? atol_ : atol_fsa_);
}

void Solver::set_absolute_tolerance_fsa(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    atol_fsa_ = static_cast<realtype>(atol);

    if (get_init_done()) {
        apply_sensitivity_tolerances();
    }
}

double Solver::get_relative_tolerance_b() const {
    return static_cast<double>(isnan(rtolB_) ? rtol_ : rtolB_);
}

void Solver::set_relative_tolerance_b(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtolB_ = static_cast<realtype>(rtol);

    if (get_init_done()) {
        apply_sensitivity_tolerances();
    }
}

double Solver::get_absolute_tolerance_b() const {
    return static_cast<double>(isnan(atolB_) ? atol_ : atolB_);
}

void Solver::set_absolute_tolerance_b(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    atolB_ = static_cast<realtype>(atol);

    if (get_init_done()) {
        apply_sensitivity_tolerances();
    }
}

double Solver::get_relative_tolerance_quadratures() const {
    return static_cast<double>(quad_rtol_);
}

void Solver::set_relative_tolerance_quadratures(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    quad_rtol_ = static_cast<realtype>(rtol);

    if (sensi_meth_ != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int)solver_memory_B_.size(); ++iMem)
        if (solver_memory_B_.at(iMem))
            apply_quad_tolerances_asa(iMem);
}

double Solver::get_absolute_tolerance_quadratures() const {
    return static_cast<double>(quad_atol_);
}

void Solver::set_absolute_tolerance_quadratures(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    quad_atol_ = static_cast<realtype>(atol);

    if (sensi_meth_ != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int)solver_memory_B_.size(); ++iMem)
        if (solver_memory_B_.at(iMem))
            apply_tolerances_asa(iMem);
}

double Solver::get_steady_state_tolerance_factor() const {
    return static_cast<double>(ss_tol_factor_);
}

void Solver::set_steady_state_tolerance_factor(double const factor) {
    if (factor < 0)
        throw AmiException("ss_tol_factor must be a non-negative number");

    ss_tol_factor_ = static_cast<realtype>(factor);
}

double Solver::get_relative_tolerance_steady_state() const {
    return static_cast<double>(
        isnan(ss_rtol_) ? rtol_ * ss_tol_factor_ : ss_rtol_
    );
}

void Solver::set_relative_tolerance_steady_state(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    ss_rtol_ = static_cast<realtype>(rtol);
}

double Solver::get_absolute_tolerance_steady_state() const {
    return static_cast<double>(
        isnan(ss_atol_) ? atol_ * ss_tol_factor_ : ss_atol_
    );
}

void Solver::set_absolute_tolerance_steady_state(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    ss_atol_ = static_cast<realtype>(atol);
}

double Solver::get_steady_state_sensi_tolerance_factor() const {
    return static_cast<double>(ss_tol_sensi_factor_);
}

void Solver::set_steady_state_sensi_tolerance_factor(
    double const ss_tol_sensi_factor
) {
    if (ss_tol_sensi_factor < 0)
        throw AmiException("ss_tol_sensi_factor must be a non-negative number");

    ss_tol_sensi_factor_ = static_cast<realtype>(ss_tol_sensi_factor);
}

double Solver::get_relative_tolerance_steady_state_sensi() const {
    return static_cast<double>(
        isnan(ss_rtol_sensi_) ? rtol_ * ss_tol_sensi_factor_ : ss_rtol_sensi_
    );
}

void Solver::set_relative_tolerance_steady_state_sensi(double const rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    ss_rtol_sensi_ = static_cast<realtype>(rtol);
}

double Solver::get_absolute_tolerance_steady_state_sensi() const {
    return static_cast<double>(
        isnan(ss_atol_sensi_) ? atol_ * ss_tol_sensi_factor_ : ss_atol_sensi_
    );
}

void Solver::set_absolute_tolerance_steady_state_sensi(double const atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    ss_atol_sensi_ = static_cast<realtype>(atol);
}

long int Solver::get_max_steps() const { return maxsteps_; }

double Solver::get_max_time() const { return maxtime_.count(); }

void Solver::set_max_time(double const maxtime) {
    maxtime_ = std::chrono::duration<double>(maxtime);
}

void Solver::start_timer() const { simulation_timer_.reset(); }

bool Solver::time_exceeded(int const interval) const {
    thread_local int eval_counter = 0;

    // 0 means infinite time
    if (maxtime_.count() == 0)
        return false;

    if (++eval_counter % interval)
        return false;

    eval_counter = 0;
    auto elapsed_s = simulation_timer_.elapsed_seconds();
    return std::chrono::duration<double>(elapsed_s) > maxtime_;
}

void Solver::set_max_steps(long int const maxsteps) {
    if (maxsteps <= 0)
        throw AmiException("maxsteps must be a positive number");

    maxsteps_ = maxsteps;
    if (get_adj_init_done())
        reset_mutable_memory(nx(), nplist(), nquad());
}

long int Solver::get_max_steps_backward_problem() const { return maxstepsB_; }

void Solver::set_max_steps_backward_problem(long int const maxsteps) {
    if (maxsteps < 0)
        throw AmiException("maxsteps must be a non-negative number");

    maxstepsB_ = maxsteps;
}

LinearMultistepMethod Solver::get_linear_multistep_method() const {
    return lmm_;
}

void Solver::set_linear_multistep_method(LinearMultistepMethod const lmm) {
    if (solver_memory_)
        reset_mutable_memory(nx(), nplist(), nquad());
    lmm_ = lmm;
}

NonlinearSolverIteration Solver::get_non_linear_solver_iteration() const {
    return iter_;
}

void Solver::set_non_linear_solver_iteration(
    NonlinearSolverIteration const iter
) {
    if (solver_memory_)
        reset_mutable_memory(nx(), nplist(), nquad());
    iter_ = iter;
}

InterpolationType Solver::get_interpolation_type() const {
    return interp_type_;
}

void Solver::set_interpolation_type(InterpolationType const interpType) {
    if (!solver_memory_B_.empty())
        reset_mutable_memory(nx(), nplist(), nquad());
    interp_type_ = interpType;
}

int Solver::get_state_ordering() const { return ordering_; }

void Solver::set_state_ordering(int ordering) {
    ordering_ = ordering;
    if (solver_memory_ && linsol_ == LinearSolver::KLU) {
        auto klu = dynamic_cast<SUNLinSolKLU*>(linear_solver_.get());
        klu->set_ordering(static_cast<SUNLinSolKLU::StateOrdering>(ordering));
        klu = dynamic_cast<SUNLinSolKLU*>(linear_solver_B_.get());
        klu->set_ordering(static_cast<SUNLinSolKLU::StateOrdering>(ordering));
    }
#ifdef SUNDIALS_SUPERLUMT
    if (solverMemory && linsol == LinearSolver::SuperLUMT) {
        auto klu = dynamic_cast<SUNLinSolSuperLUMT*>(linearSolver.get());
        klu->setOrdering(
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(ordering)
        );
        klu = dynamic_cast<SUNLinSolSuperLUMT*>(linearSolverB.get());
        klu->setOrdering(
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(ordering)
        );
    }
#endif
}

bool Solver::get_stability_limit_flag() const { return stldet_; }

void Solver::set_stability_limit_flag(bool const stldet) {
    stldet_ = stldet;
    if (solver_memory_) {
        set_stab_lim_det(stldet);
        for (int iMem = 0; iMem < (int)solver_memory_B_.size(); ++iMem)
            if (solver_memory_B_.at(iMem))
                set_stab_lim_det_b(iMem, stldet);
    }
}

LinearSolver Solver::get_linear_solver() const { return linsol_; }

void Solver::set_linear_solver(LinearSolver const linsol) {
    if (solver_memory_)
        reset_mutable_memory(nx(), nplist(), nquad());
    linsol_ = linsol;
}

InternalSensitivityMethod Solver::get_internal_sensitivity_method() const {
    return ism_;
}

void Solver::set_internal_sensitivity_method(
    InternalSensitivityMethod const ism
) {
    if (solver_memory_)
        reset_mutable_memory(nx(), nplist(), nquad());
    ism_ = ism;
}

RDataReporting Solver::get_return_data_reporting_mode() const {
    return rdata_mode_;
}

void Solver::set_return_data_reporting_mode(RDataReporting const rdrm) {
    if (rdrm == RDataReporting::residuals
        && sensi_meth_ == SensitivityMethod::adjoint)
        throw AmiException(
            "Adjoint Sensitivity Analysis cannot report "
            "residuals!"
        );
    rdata_mode_ = rdrm;
}

void Solver::initialize_non_linear_solver_sens(Model const* model) const {
    switch (iter_) {
    case NonlinearSolverIteration::newton:
        switch (ism_) {
        case InternalSensitivityMethod::staggered:
        case InternalSensitivityMethod::simultaneous:
            non_linear_solver_sens_ = std::make_unique<SUNNonLinSolNewton>(
                1 + model->nplist(), x_.get_nvector()
            );
            break;
        case InternalSensitivityMethod::staggered1:
            non_linear_solver_sens_
                = std::make_unique<SUNNonLinSolNewton>(x_.get_nvector());
            break;
        default:
            throw AmiException(
                "Unsupported internal sensitivity method selected: %d", ism_
            );
        }
        break;
    case NonlinearSolverIteration::fixedpoint:
        switch (ism_) {
        case InternalSensitivityMethod::staggered:
        case InternalSensitivityMethod::simultaneous:
            non_linear_solver_sens_ = std::make_unique<SUNNonLinSolFixedPoint>(
                1 + model->nplist(), x_.get_nvector()
            );
            break;
        case InternalSensitivityMethod::staggered1:
            non_linear_solver_sens_
                = std::make_unique<SUNNonLinSolFixedPoint>(x_.get_nvector());
            break;
        default:
            throw AmiException(
                "Unsupported internal sensitivity method selected: %d", ism_
            );
        }
        break;
    default:
        throw AmiException(
            "Invalid non-linear solver specified (%d).", static_cast<int>(iter_)
        );
    }

    set_non_linear_solver_sens();
}

void Solver::set_err_handler_fn() const {
    SUNContext_ClearErrHandlers(sunctx_);
    auto sunerr = SUNContext_PushErrHandler(
        sunctx_, wrapErrHandlerFn,
        reinterpret_cast<void*>(const_cast<Solver*>(this))
    );
    if (sunerr)
        throw AmiException(
            "Error setting error handler function: %s", SUNGetErrMsg(sunerr)
        );
}

int Solver::nplist() const { return sx_.size(); }

int Solver::nx() const { return x_.size(); }

int Solver::nquad() const { return xQB_.size(); }

bool Solver::get_init_done() const { return initialized_; }

bool Solver::get_sens_init_done() const { return sens_initialized_; }

bool Solver::get_adj_init_done() const { return adj_initialized_; }

bool Solver::get_init_done_b(int const which) const {
    return static_cast<int>(initializedB_.size()) > which
           && initializedB_.at(which);
}

bool Solver::get_quad_init_done_b(int const which) const {
    return static_cast<int>(initializedQB_.size()) > which
           && initializedQB_.at(which);
}

bool Solver::get_quad_init_done() const { return quad_initialized_; }

void Solver::set_init_done() const { initialized_ = true; }

void Solver::set_sens_init_done() const { sens_initialized_ = true; }

void Solver::set_adj_init_done() const { adj_initialized_ = true; }

void Solver::set_init_done_b(int const which) const {
    if (which >= static_cast<int>(initializedB_.size()))
        initializedB_.resize(which + 1, false);
    initializedB_.at(which) = true;
}

void Solver::set_quad_init_done_b(int const which) const {
    if (which >= static_cast<int>(initializedQB_.size()))
        initializedQB_.resize(which + 1, false);
    initializedQB_.at(which) = true;
}

void Solver::set_quad_init_done() const { quad_initialized_ = true; }

void Solver::switch_forward_sensis_off() const { sens_toggle_off(); }

realtype Solver::get_cpu_time() const { return cpu_time_; }

realtype Solver::get_cpu_time_b() const { return cpu_time_b_; }

void Solver::reset_mutable_memory(
    int const nx, int const nplist, int const nquad
) const {
    solver_memory_ = nullptr;
    initialized_ = false;
    adj_initialized_ = false;
    sens_initialized_ = false;
    quad_initialized_ = false;
    solver_was_called_F_ = false;
    solver_was_called_B_ = false;

    x_ = AmiVector(nx, sunctx_);
    dx_ = AmiVector(nx, sunctx_);
    sx_ = AmiVectorArray(nx, nplist, sunctx_);
    sdx_ = AmiVectorArray(nx, nplist, sunctx_);

    dky_ = AmiVector(nx, sunctx_);

    xB_ = AmiVector(nx, sunctx_);
    dxB_ = AmiVector(nx, sunctx_);
    xQB_ = AmiVector(nquad, sunctx_);
    xQ_ = AmiVector(nx, sunctx_);

    solver_memory_B_.clear();
    initializedB_.clear();
    initializedQB_.clear();
}

void Solver::write_solution(
    realtype& t, AmiVector& x, AmiVector& dx, AmiVectorArray& sx, AmiVector& xQ
) const {
    t = get_t();
    if (quad_initialized_)
        xQ.copy(get_quadrature(t));
    if (sens_initialized_)
        sx.copy(get_state_sensitivity(t));
    x.copy(get_state(t));
    dx.copy(get_derivative_state(t));
}

void Solver::write_solution(
    realtype& t, AmiVector& x, AmiVector& dx, AmiVectorArray& sx
) const {
    t = get_t();
    if (sens_initialized_)
        sx.copy(get_state_sensitivity(t));
    x.copy(get_state(t));
    dx.copy(get_derivative_state(t));
}

void Solver::write_solution(SolutionState& sol) const {
    sol.t = get_t();
    if (sens_initialized_)
        sol.sx.copy(get_state_sensitivity(sol.t));
    sol.x.copy(get_state(sol.t));
    sol.dx.copy(get_derivative_state(sol.t));
}

void Solver::write_solution(realtype const t, SolutionState& sol) const {
    sol.t = t;
    if (sens_initialized_)
        sol.sx.copy(get_state_sensitivity(sol.t));
    sol.x.copy(get_state(sol.t));
    sol.dx.copy(get_derivative_state(sol.t));
}

void Solver::write_solution_b(
    realtype& t, AmiVector& xB, AmiVector& dxB, AmiVector& xQB, int const which
) const {
    t = get_t();
    xB.copy(get_adjoint_state(which, t));
    dxB.copy(get_adjoint_derivative_state(which, t));
    xQB.copy(get_adjoint_quadrature(which, t));
}

AmiVector const& Solver::get_state(realtype const t) const {
    if (t == t_)
        return x_;

    if (solver_was_called_F_)
        get_dky(t, 0);

    return dky_;
}

AmiVector const& Solver::get_derivative_state(realtype const t) const {
    if (t == t_)
        return dx_;

    if (solver_was_called_F_)
        get_dky(t, 1);

    return dky_;
}

AmiVectorArray const& Solver::get_state_sensitivity(realtype const t) const {
    if (sens_initialized_ && solver_was_called_F_) {
        if (t == t_) {
            get_sens();
        } else {
            get_sens_dky(t, 0);
        }
    }
    return sx_;
}

AmiVector const&
Solver::get_adjoint_state(int const which, realtype const t) const {
    if (adj_initialized_) {
        if (solver_was_called_B_) {
            if (t == t_) {
                get_b(which);
                return xB_;
            }
            get_dky_b(t, 0, which);
        }
    } else {
        dky_.zero();
    }
    return dky_;
}

AmiVector const&
Solver::get_adjoint_derivative_state(int const which, realtype const t) const {
    if (adj_initialized_) {
        if (solver_was_called_B_) {
            if (t == t_) {
                get_b(which);
                return dxB_;
            }
            get_dky_b(t, 1, which);
        }
    } else {
        dky_.zero();
    }
    return dky_;
}

AmiVector const&
Solver::get_adjoint_quadrature(int const which, realtype const t) const {
    if (adj_initialized_) {
        if (solver_was_called_B_) {
            if (t == t_) {
                get_quad_b(which);
                return xQB_;
            }
            get_quad_dky_b(t, 0, which);
        }
    } else {
        xQB_.zero();
    }
    return xQB_;
}

AmiVector const& Solver::get_quadrature(realtype t) const {
    if (quad_initialized_) {
        if (solver_was_called_F_) {
            if (t == t_) {
                get_quad(t);
                return xQ_;
            }
            get_quad_dky(t, 0);
        }
    } else {
        xQ_.zero();
    }
    return xQ_;
}

realtype Solver::get_t() const { return t_; }
} // namespace amici
