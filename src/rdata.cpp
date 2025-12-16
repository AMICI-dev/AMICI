#include "amici/rdata.h"

#include "amici/backwardproblem.h"
#include "amici/edata.h"
#include "amici/forwardproblem.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/steadystateproblem.h"
#include "amici/symbolic_functions.h"
#include "amici/vector.h"

#include <algorithm>
#include <cmath>

namespace amici {
using std::isnan;

ReturnData::ReturnData(Solver const& solver, Model const& model)
    : ReturnData(
          model.get_timepoints(),
          ModelDimensions(static_cast<ModelDimensions const&>(model)),
          model.n_max_event(), solver.get_newton_max_steps(),
          model.get_parameter_list(), model.get_parameter_scale(),
          model.get_second_order_mode(), solver.get_sensitivity_order(),
          solver.get_sensitivity_method(),
          solver.get_return_data_reporting_mode(), model.has_quadratic_llh(),
          model.get_add_sigma_residuals(), model.get_minimum_sigma_residuals(),
          model.get_free_parameter_ids(), model.get_free_parameter_names(),
          model.get_fixed_parameter_ids(), model.get_fixed_parameter_names(),
          model.get_state_ids(), model.get_state_names(),
          model.get_state_ids_solver(), model.get_state_names_solver(),
          model.get_observable_ids(), model.get_observable_names(),
          model.get_expression_ids(), model.get_expression_names()
      ) {}

ReturnData::ReturnData(
    std::vector<realtype> ts_, ModelDimensions const& model_dimensions_,
    int nmaxevent_, int newton_maxsteps_, std::vector<int> plist_,
    std::vector<ParameterScaling> pscale_, SecondOrderMode o2mode_,
    SensitivityOrder sensi_, SensitivityMethod sensi_meth_,
    RDataReporting rdrm_, bool quadratic_llh_, bool sigma_res_,
    realtype sigma_offset_,
    std::span<std::string_view const> free_parameter_ids_,
    std::span<std::string_view const> free_parameter_names_,
    std::span<std::string_view const> fixed_parameter_ids_,
    std::span<std::string_view const> fixed_parameter_names_,
    std::span<std::string_view const> state_ids_,
    std::span<std::string_view const> state_names_,
    std::span<std::string_view const> state_ids_solver_,
    std::span<std::string_view const> state_names_solver_,
    std::span<std::string_view const> observable_ids_,
    std::span<std::string_view const> observable_names_,
    std::span<std::string_view const> expression_ids_,
    std::span<std::string_view const> expression_names_
)
    : ModelDimensions(model_dimensions_)
    , ts(std::move(ts_))
    , nplist(plist_.size())
    , nmaxevent(nmaxevent_)
    , nt(ts.size())
    , newton_maxsteps(newton_maxsteps_)
    , pscale(std::move(pscale_))
    , o2mode(o2mode_)
    , sensi(sensi_)
    , sensi_meth(sensi_meth_)
    , rdata_reporting(rdrm_)
    , sigma_res(sigma_res_)
    , plist(plist_)
    , free_parameter_ids(free_parameter_ids_)
    , free_parameter_names(free_parameter_names_)
    , fixed_parameter_ids(fixed_parameter_ids_)
    , fixed_parameter_names(fixed_parameter_names_)
    , state_ids(state_ids_)
    , state_names(state_names_)
    , state_ids_solver(state_ids_solver_)
    , state_names_solver(state_names_solver_)
    , observable_ids(observable_ids_)
    , observable_names(observable_names_)
    , expression_ids(expression_ids_)
    , expression_names(expression_names_)
    , sigma_offset(sigma_offset_)
    , nroots_(ne) {
    model_dimensions_.validate();

    switch (rdata_reporting) {
    case RDataReporting::full:
        initialize_full_reporting(quadratic_llh_);
        break;

    case RDataReporting::residuals:
        initialize_residual_reporting(quadratic_llh_);
        break;

    case RDataReporting::likelihood:
        initialize_likelihood_reporting(quadratic_llh_);
        break;

    case RDataReporting::observables_likelihood:
        initialize_observables_likelihood_reporting(quadratic_llh_);
        break;
    }
}

void ReturnData::initialize_likelihood_reporting(bool enable_fim) {
    llh = get_nan();
    chi2 = get_nan();
    if (sensi >= SensitivityOrder::first
        && sensi_meth != SensitivityMethod::none) {
        sllh.resize(nplist, get_nan());
        if (sensi >= SensitivityOrder::second)
            s2llh.resize(nplist * (nJ - 1), get_nan());

        if ((sensi_meth == SensitivityMethod::forward
             || sensi >= SensitivityOrder::second)
            && enable_fim)
            FIM.resize(nplist * nplist, 0.0);
    }
}

void ReturnData::initialize_observables_likelihood_reporting(bool enable_fim) {
    initialize_likelihood_reporting(enable_fim);

    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);

    if ((sensi_meth == SensitivityMethod::forward
         && sensi >= SensitivityOrder::first)
        || sensi >= SensitivityOrder::second) {

        sy.resize(nt * ny * nplist, 0.0);
        ssigmay.resize(nt * ny * nplist, 0.0);
    }
}

void ReturnData::initialize_residual_reporting(bool enable_res) {
    y.resize(nt * ny, 0.0);
    sigmay.resize(nt * ny, 0.0);
    if (enable_res)
        res.resize((sigma_res ? 2 : 1) * nt * nytrue, 0.0);

    if ((sensi_meth == SensitivityMethod::forward
         && sensi >= SensitivityOrder::first)
        || sensi >= SensitivityOrder::second) {

        sy.resize(nt * ny * nplist, 0.0);
        ssigmay.resize(nt * ny * nplist, 0.0);
        if (enable_res)
            sres.resize((sigma_res ? 2 : 1) * nt * nytrue * nplist, 0.0);
    }
}

void ReturnData::initialize_full_reporting(bool quadratic_llh) {

    initialize_likelihood_reporting(quadratic_llh);
    initialize_residual_reporting(quadratic_llh);

    xdot.resize(nx_solver, get_nan());

    J.resize(nx_solver * nx_solver, get_nan());

    // initialize with 0.0, so we only need to write non-zero values
    z.resize(nmaxevent * nz, 0.0);
    sigmaz.resize(nmaxevent * nz, 0.0);

    rz.resize(nmaxevent * nz, 0.0);
    x.resize(nt * nx_rdata, 0.0);
    w.resize(nt * nw, 0.0);

    preeq_numsteps.resize(3, 0);
    preeq_status.resize(3, SteadyStateStatus::not_run);
    posteq_numsteps.resize(3, 0);
    posteq_status.resize(3, SteadyStateStatus::not_run);

    if (nt > 0) {
        numsteps.resize(nt, 0);
        num_rhs_evals.resize(nt, 0);
        num_err_test_fails.resize(nt, 0);
        num_non_lin_solv_conv_fails.resize(nt, 0);
        order.resize(nt, 0);

        if (sensi_meth == SensitivityMethod::adjoint
            && sensi >= SensitivityOrder::first) {
            numsteps_b.resize(nt, 0);
            num_rhs_evals_b.resize(nt, 0);
            num_err_test_fails_b.resize(nt, 0);
            num_non_lin_solv_conv_fails_b.resize(nt, 0);
        }
    }

    x0.resize(nx_rdata, get_nan());
    x_ss.resize(nx_rdata, get_nan());

    if (sensi >= SensitivityOrder::first) {
        sx0.resize(nx_rdata * nplist, get_nan());
        sx_ss.resize(nx_rdata * nplist, get_nan());

        if (sensi_meth == SensitivityMethod::forward
            || sensi >= SensitivityOrder::second) {
            // for second order we can fill in from the augmented states
            sx.resize(nt * nx_rdata * nplist, 0.0);
            sz.resize(nmaxevent * nz * nplist, 0.0);
            srz.resize(nmaxevent * nz * nplist, 0.0);
        }

        ssigmaz.resize(nmaxevent * nz * nplist, 0.0);
        if (sensi >= SensitivityOrder::second
            && sensi_meth == SensitivityMethod::forward)
            s2rz.resize(nmaxevent * nztrue * nplist * nplist, 0.0);
    }
}

void ReturnData::process_simulation_objects(
    ForwardProblem const* fwd, BackwardProblem const* bwd, Model& model,
    Solver const& solver, ExpData const* edata
) {
    ModelContext mc(&model);

    process_solver(solver);

    SteadyStateProblem const* preeq = nullptr;
    SteadyStateProblem const* posteq = nullptr;
    SteadyStateBackwardProblem const* preeq_bwd = nullptr;
    SteadyStateBackwardProblem const* posteq_bwd = nullptr;
    if (fwd) {
        preeq = fwd->get_preeq_problem();
        posteq = fwd->get_posteq_problem();
        if (bwd) {
            preeq_bwd = bwd->get_preeq_bwd_problem();
            posteq_bwd = bwd->get_posteq_bwd_problem();
        }
    }
    if (preeq)
        process_pre_equilibration(*preeq, preeq_bwd, model);

    if (fwd)
        process_forward_problem(*fwd, model, edata);
    else
        invalidate(0);

    if (posteq)
        process_post_equilibration(*posteq, posteq_bwd, model, edata);

    if (edata && !edata->fixed_parameters_pre_equilibration.empty() && fwd
        && !fwd->was_preequilibrated() && preeq) {
        // failure during preequilibration
        store_jacobian_and_derivative(*preeq, model);
    } else if (fwd && !posteq)
        store_jacobian_and_derivative(*fwd, model);
    else if (posteq)
        store_jacobian_and_derivative(*posteq, model);

    if (fwd && bwd)
        process_backward_problem(*fwd, *bwd, preeq, preeq_bwd, model);
    else if (solver.computing_asa())
        invalidate_sllh();

    apply_chain_rule_factor_to_simulation_results(model);
}

void ReturnData::process_pre_equilibration(
    SteadyStateProblem const& preeq,
    SteadyStateBackwardProblem const* preeq_bwd, Model& model
) {
    auto const simulation_state = preeq.get_final_simulation_state();
    model.set_model_state(simulation_state.mod);

    if (!x_ss.empty()) {
        model.fx_rdata(x_ss, simulation_state.sol.x);
    }
    if (!sx_ss.empty() && sensi >= SensitivityOrder::first) {
        model.fsx_rdata(sx_ss, simulation_state.sol.sx, simulation_state.sol.x);
    }
    /* Get cpu time for pre-equilibration in milliseconds */
    preeq_cpu_time = preeq.get_cpu_time();
    if (preeq_bwd) {
        preeq_cpu_time_b = preeq_bwd->get_cpu_time_b();
        preeq_numsteps_b = preeq_bwd->get_num_steps_b();
    }
    preeq_wrms = preeq.get_residual_norm();
    preeq_status = preeq.get_steady_state_status();
    if (preeq_status[1] == SteadyStateStatus::success)
        preeq_t = preeq.get_steady_state_time();
    if (!preeq_numsteps.empty())
        write_slice(preeq.get_num_steps(), preeq_numsteps);
}

void ReturnData::process_post_equilibration(
    SteadyStateProblem const& posteq,
    SteadyStateBackwardProblem const* posteq_bwd, Model& model,
    ExpData const* edata
) {
    for (int it = 0; it < nt; it++) {
        auto const t = model.get_timepoint(it);
        if (std::isinf(t)) {
            auto const& simulation_state = posteq.get_final_simulation_state();
            model.set_model_state(simulation_state.mod);
            get_data_output(it, model, simulation_state.sol, edata);
        }
    }
    /* Get cpu time for Newton solve in milliseconds */
    posteq_cpu_time = posteq.get_cpu_time();
    if (posteq_bwd) {
        posteq_cpu_time_b = posteq_bwd->get_cpu_time_b();
        posteq_numsteps_b = posteq_bwd->get_num_steps_b();
    }
    posteq_wrms = posteq.get_residual_norm();
    posteq_status = posteq.get_steady_state_status();
    if (posteq_status[1] == SteadyStateStatus::success)
        posteq_t = posteq.get_steady_state_time();
    if (!posteq_numsteps.empty())
        write_slice(posteq.get_num_steps(), posteq_numsteps);
}

void ReturnData::process_forward_problem(
    ForwardProblem const& fwd, Model& model, ExpData const* edata
) {
    if (edata)
        initialize_objective_function(model.has_quadratic_llh());

    auto const& initialState = fwd.get_initial_simulation_state();
    if (initialState.sol.x.size() == 0 && model.nx_solver > 0)
        return; // if x wasn't set forward problem failed during initialization

    model.set_model_state(initialState.mod);

    if (!x0.empty()) {
        model.fx_rdata(x0, initialState.sol.x);
    }

    if (!sx0.empty()) {
        model.fsx_rdata(sx0, initialState.sol.sx, initialState.sol.x);
    }

    // process timepoint data
    realtype tf = fwd.get_final_time();
    for (int it = 0; it < model.nt(); it++) {
        if (model.get_timepoint(it) <= tf) {
            auto const simulation_state
                = fwd.get_simulation_state_timepoint(it);
            model.set_model_state(simulation_state.mod);
            get_data_output(it, model, simulation_state.sol, edata);
        } else {
            // check for integration failure but consider postequilibration
            if (!std::isinf(model.get_timepoint(it)))
                invalidate(it);
        }
    }

    // process event data
    if (nz > 0) {
        auto const& discontinuities = fwd.get_discontinuities();
        Expects(
            static_cast<int>(discontinuities.size())
            == fwd.get_event_counter() + 1
        );
        for (int iroot = 0; iroot <= fwd.get_event_counter(); iroot++) {
            auto const simulation_state = fwd.get_simulation_state_event(iroot);
            model.set_model_state(simulation_state.mod);
            get_event_output(
                discontinuities.at(iroot).root_info, model,
                simulation_state.sol, edata
            );
        }
    }
}

void ReturnData::get_data_output(
    int it, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (!x.empty()) {
        model.fx_rdata(slice(x, it, nx_rdata), sol.x);
    }
    if (!w.empty())
        model.get_expression(slice(w, it, nw), ts[it], sol.x);
    if (!y.empty())
        model.get_observable(slice(y, it, ny), ts[it], sol.x);
    if (!sigmay.empty())
        model.get_observable_sigma(slice(sigmay, it, ny), it, edata);

    if (edata) {
        if (!isnan(llh))
            model.add_observable_objective(llh, it, sol.x, *edata);
        fres(it, model, sol, *edata);
        fchi2(it, *edata);
    }

    if (sensi >= SensitivityOrder::first && nplist > 0) {

        if (sensi_meth == SensitivityMethod::forward) {
            get_data_sensis_fsa(it, model, sol, edata);
        } else if (edata && !sllh.empty()) {
            model.add_partial_observable_objective_sensitivity(
                sllh, s2llh, it, sol.x, *edata
            );
        }

        if (!ssigmay.empty())
            model.get_observable_sigma_sensitivity(
                slice(ssigmay, it, nplist * ny), slice(sy, it, nplist * ny), it,
                edata
            );
    }
}

void ReturnData::get_data_sensis_fsa(
    int it, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (!sx.empty()) {
        model.fsx_rdata(
            gsl::span<realtype>(
                &sx.at(it * nplist * nx_rdata), nplist * nx_rdata
            ),
            sol.sx, sol.x
        );
    }

    if (!sy.empty()) {
        model.get_observable_sensitivity(
            slice(sy, it, nplist * ny), ts[it], sol.x, sol.sx
        );
    }

    if (edata) {
        if (!sllh.empty())
            model.add_observable_objective_sensitivity(
                sllh, s2llh, it, sol.x, sol.sx, *edata
            );
        fsres(it, model, sol, *edata);
        fFIM(it, model, sol, *edata);
    }
}

void ReturnData::get_event_output(
    std::vector<int> const& rootidx, Model& model, SolutionState const& sol,
    ExpData const* edata
) {

    for (int ie = 0; ie < ne; ie++) {
        if (rootidx.at(ie) != 1 || nroots_.at(ie) >= nmaxevent)
            continue;

        /* get event output */
        if (!z.empty())
            model.get_event(slice(z, nroots_.at(ie), nz), ie, sol.t, sol.x);
        /* if called from fillEvent at last timepoint,
         then also get the root function value */
        if (sol.t == model.get_timepoint(nt - 1))
            if (!rz.empty())
                model.get_event_regularization(
                    slice(rz, nroots_.at(ie), nz), ie, sol.t, sol.x
                );

        if (edata) {
            if (!sigmaz.empty())
                model.get_event_sigma(
                    slice(sigmaz, nroots_.at(ie), nz), ie, nroots_.at(ie),
                    sol.t, edata
                );
            if (!isnan(llh))
                model.add_event_objective(
                    llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );

            /* if called from fillEvent at last timepoint,
               add regularization based on rz */
            if (sol.t == model.get_timepoint(nt - 1) && !isnan(llh)) {
                model.add_event_objective_regularization(
                    llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );
            }
        }

        if (sensi >= SensitivityOrder::first) {
            if (sensi_meth == SensitivityMethod::forward) {
                get_event_sensis_fsa(ie, model, sol, edata);
            } else if (edata && !sllh.empty()) {
                model.add_partial_event_objective_sensitivity(
                    sllh, s2llh, ie, nroots_.at(ie), sol.t, sol.x, *edata
                );
            }
        }
        nroots_.at(ie)++;
    }
}

void ReturnData::get_event_sensis_fsa(
    int ie, Model& model, SolutionState const& sol, ExpData const* edata
) {
    if (sol.t == model.get_timepoint(nt - 1)) {
        // call from fillEvent at last timepoint
        if (!sz.empty())
            model.get_unobserved_event_sensitivity(
                slice(sz, nroots_.at(ie), nz * nplist), ie
            );
        if (!srz.empty())
            model.get_event_regularization_sensitivity(
                slice(srz, nroots_.at(ie), nz * nplist), ie, sol.t, sol.x,
                sol.sx
            );
    } else if (!sz.empty()) {
        model.get_event_sensitivity(
            slice(sz, nroots_.at(ie), nz * nplist), ie, sol.t, sol.x, sol.sx
        );
    }

    if (edata && !sllh.empty()) {
        model.add_event_objective_sensitivity(
            sllh, s2llh, ie, nroots_.at(ie), sol.t, sol.x, sol.sx, *edata
        );
    }
}

void ReturnData::process_backward_problem(
    ForwardProblem const& fwd, BackwardProblem const& bwd,
    SteadyStateProblem const* preeq,
    SteadyStateBackwardProblem const* preeq_bwd, Model& model
) {
    if (sllh.empty())
        return;

    auto simulation_state = fwd.get_initial_simulation_state();
    model.set_model_state(simulation_state.mod);

    std::vector<realtype> llhS0(model.nJ * model.nplist(), 0.0);
    // Adjoint quadratures before pre-equilibration
    auto xQB = bwd.get_adjoint_quadrature_pre_preeq();

    if (preeq_bwd && preeq_bwd->has_quadrature()) {
        // Pre-equilibration with ASA steady-state shortcuts

        // If pre-equilibration is run in adjoint mode, the scalar product of
        // sx0 with its adjoint counterpart (see handleSx0Forward()) is not
        // necessary:
        // the actual simulation is "extended" by the pre-equilibration time.
        // At initialization (at t=-inf), the adjoint state is in steady state
        // (= 0) and so is the scalar product.
        // Instead of the scalar product, the quadratures xQB from
        // pre-equilibration contribute to the gradient
        // (see example notebook on equilibration for further documentation).

        // Adjoint quadratures after pre-equilibration
        auto const& xQB_post_preeq = preeq_bwd->get_adjoint_quadrature();
        for (int ip = 0; ip < model.nplist(); ++ip)
            xQB[ip] += xQB_post_preeq.at(ip);

        handle_sx0_backward(
            model, preeq->get_state_sensitivity(), bwd.get_adjoint_state(),
            llhS0
        );
    } else if (preeq
               && preeq->get_steady_state_status()[1]
                      != SteadyStateStatus::not_run) {
        // Pre-equilibration with ASA backward simulation

        // Adjoint quadratures after pre-equilibration
        xQB = bwd.get_adjoint_quadrature();

        handle_sx0_backward(
            model, preeq->get_state_sensitivity(), bwd.get_adjoint_state(),
            llhS0
        );

    } else {
        // No pre-equilibration, or pre-equilibration with FSA
        handle_sx0_forward(
            model, simulation_state.sol, llhS0,
            bwd.get_adjoint_state_pre_preeq()
        );
    }

    for (int iJ = 0; iJ < model.nJ; iJ++) {
        for (int ip = 0; ip < model.nplist(); ip++) {
            if (iJ == 0) {
                sllh.at(ip) -= llhS0[ip] + xQB[ip * model.nJ];
            } else {
                s2llh.at(iJ - 1 + ip * (model.nJ - 1))
                    -= llhS0[ip + iJ * model.nplist()]
                       + xQB[iJ + ip * model.nJ];
            }
        }
    }
}

void ReturnData::handle_sx0_backward(
    Model const& model, AmiVectorArray const& sx0, AmiVector const& xB,
    std::vector<realtype>& llhS0
) const {

    // Add the contribution for sx0 from preequilibration. If backward
    // preequilibration was done by simulation due to a singular Jacobian,
    // xB is not necessarily 0 and we may get a non-zero contribution here.
    for (int ip = 0; ip < model.nplist(); ++ip) {
        llhS0[ip] = 0.0;
        for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
            llhS0[ip] += xB.at(ix) * sx0.at(ix, ip);
        }
    }
}

void ReturnData::handle_sx0_forward(
    Model const& model, SolutionState const& sol, std::vector<realtype>& llhS0,
    AmiVector const& xB
) const {
    /* If preequilibration is run in forward mode or is not needed, then adjoint
       sensitivity analysis still needs the state sensitivities at t=0 (sx0),
       to compute the gradient. For each parameter, the scalar product of sx0
       with its adjoint counterpart contributes to the gradient. */
    for (int iJ = 0; iJ < model.nJ; iJ++) {
        if (iJ == 0) {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip] += xB[ix] * sol.sx.at(ix, ip);
                }
            }
        } else {
            for (int ip = 0; ip < model.nplist(); ++ip) {
                llhS0[ip + iJ * model.nplist()] = 0.0;
                for (int ix = 0; ix < model.nxtrue_solver; ++ix) {
                    llhS0[ip + iJ * model.nplist()]
                        += xB[ix + iJ * model.nxtrue_solver] * sol.sx.at(ix, ip)
                           + xB[ix]
                                 * sol.sx.at(ix + iJ * model.nxtrue_solver, ip);
                }
            }
        }
    }
}

void ReturnData::process_solver(Solver const& solver) {
    using std::ranges::copy;

    cpu_time = solver.get_cpu_time();

    if (!numsteps.empty()) {
        Expects(numsteps.size() >= solver.get_num_steps().size());
        // copy instead of assignment to ensure length `nt`
        // (vector from solver may be shorter in case of integration errors)
        copy(solver.get_num_steps(), numsteps.begin());
    }

    if (!num_rhs_evals.empty()) {
        copy(solver.get_num_rhs_evals(), num_rhs_evals.begin());
    }

    if (!num_err_test_fails.empty()) {
        copy(solver.get_num_err_test_fails(), num_err_test_fails.begin());
    }

    if (!num_non_lin_solv_conv_fails.empty()) {
        copy(
            solver.get_num_non_lin_solv_conv_fails(),
            num_non_lin_solv_conv_fails.begin()
        );
    }

    if (!order.empty()) {
        copy(solver.get_last_order(), order.begin());
    }

    cpu_time_b = solver.get_cpu_time_b();

    if (!numsteps_b.empty()) {
        copy(solver.get_num_steps_b(), numsteps_b.begin());
    }

    if (!num_rhs_evals_b.empty()) {
        copy(solver.get_num_rhs_evals_b(), num_rhs_evals_b.begin());
    }

    if (!num_err_test_fails_b.empty()) {
        copy(solver.get_num_err_test_fails_b(), num_err_test_fails_b.begin());
    }

    if (!num_non_lin_solv_conv_fails_b.empty()) {
        copy(
            solver.get_num_non_lin_solv_conv_fails_b(),
            num_non_lin_solv_conv_fails_b.begin()
        );
    }
}

void ReturnData::invalidate(int const it_start) {
    if (it_start >= nt)
        return;

    invalidate_llh();
    invalidate_sllh();

    if (!x.empty())
        std::fill(x.begin() + nx_rdata * it_start, x.end(), get_nan());
    if (!y.empty())
        std::fill(y.begin() + ny * it_start, y.end(), get_nan());
    if (!w.empty())
        std::fill(w.begin() + nw * it_start, w.end(), get_nan());

    if (!sx.empty())
        std::fill(
            sx.begin() + nx_rdata * nplist * it_start, sx.end(), get_nan()
        );
    if (!sy.empty())
        std::fill(sy.begin() + ny * nplist * it_start, sy.end(), get_nan());
}

void ReturnData::invalidate_llh() {
    llh = get_nan();
    chi2 = get_nan();
}

void ReturnData::invalidate_sllh() {
    if (!sllh.empty()) {
        std::ranges::fill(sllh, get_nan());
        std::ranges::fill(s2llh, get_nan());
    }
}

void ReturnData::apply_chain_rule_factor_to_simulation_results(
    Model const& model
) {
    // chain-rule factor: multiplier for am_p
    std::vector<realtype> coefficient(nplist, 1.0);
    // Only apply `pcoefficient` to non-zero sensitivities. This allows
    // having unused parameters that are set to NAN, and still having finite
    // sensitivities. Otherwise, this would lead to NAN-sensitivities w.r.t.
    // to such parameters if they are log-scaled.
    std::vector<realtype> pcoefficient(nplist, 1.0);

    std::vector<realtype> unscaledParameters = model.get_free_parameters();
    unscale_parameters(
        unscaledParameters, model.get_parameter_scale(), unscaledParameters
    );

    std::vector<realtype> augcoefficient(np, 1.0);

    if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full) {
        for (int ip = 0; ip < np; ++ip) {
            switch (pscale[ip]) {
            case ParameterScaling::log10:
                augcoefficient.at(ip) = unscaledParameters.at(ip) * log(10);
                break;
            case ParameterScaling::ln:
                augcoefficient.at(ip) = unscaledParameters.at(ip);
                break;
            case ParameterScaling::none:
                break;
            }
        }
    }

    for (int ip = 0; ip < nplist; ++ip) {
        switch (pscale[model.plist(ip)]) {
        case ParameterScaling::log10:
            coefficient.at(ip) = log(10.0);
            pcoefficient.at(ip)
                = unscaledParameters.at(model.plist(ip)) * log(10);
            break;
        case ParameterScaling::ln:
            pcoefficient.at(ip) = unscaledParameters.at(model.plist(ip));
            break;
        case ParameterScaling::none:
            break;
        }
    }

    if (sensi >= SensitivityOrder::first) {
        // recover first order sensitivities from states for adjoint sensitivity
        // analysis
        if (sensi == SensitivityOrder::second && o2mode == SecondOrderMode::full
            && sensi_meth == SensitivityMethod::adjoint) {
            if (!sx.empty() && !x.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int ix = 0; ix < nxtrue_rdata; ++ix)
                        for (int it = 0; it < nt; ++it)
                            sx.at(ix + nxtrue_rdata * (ip + it * nplist))
                                = x.at(
                                    it * nx_rdata + nxtrue_rdata
                                    + ip * nxtrue_rdata + ix
                                );

            if (!sy.empty() && !y.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iy = 0; iy < nytrue; ++iy)
                        for (int it = 0; it < nt; ++it)
                            sy.at(iy + nytrue * (ip + it * nplist))
                                = y.at(it * ny + nytrue + ip * nytrue + iy);

            if (!sz.empty() && !z.empty())
                for (int ip = 0; ip < nplist; ++ip)
                    for (int iz = 0; iz < nztrue; ++iz)
                        for (int it = 0; it < nt; ++it)
                            sz.at(iz + nztrue * (ip + it * nplist))
                                = z.at(it * nz + nztrue + ip * nztrue + iz);
        }

        if (!sllh.empty())
            for (int ip = 0; ip < nplist; ++ip)
                if (auto&& val = sllh.at(ip); val != 0.0)
                    val *= pcoefficient.at(ip);

        if (!sres.empty())
            for (int ires = 0; ires < gsl::narrow<int>(res.size()); ++ires)
                for (int ip = 0; ip < nplist; ++ip)
                    if (auto&& val = sres.at(ires * nplist + ip); val != 0.0)
                        val *= pcoefficient.at(ip);

        if (!FIM.empty())
            for (int ip = 0; ip < nplist; ++ip)
                for (int jp = 0; jp < nplist; ++jp)
                    if (auto&& val = FIM.at(jp + ip * nplist); val != 0.0)
                        val *= pcoefficient.at(ip) * pcoefficient.at(jp);

        // apply chain rule to sensitivities
        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty()) {
                return;
            }
            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));
            for (index_type i1 = 0; i1 < gsl::narrow<index_type>(n1); ++i1)
                for (index_type ip = 0; ip < gsl::narrow<index_type>(nplist);
                     ++ip)
                    for (index_type i2 = 0; i2 < gsl::narrow<index_type>(n2);
                         ++i2)
                        if (auto&& val
                            = sens.at((i2 * nplist + ip) * stride1 + i1);
                            val != 0.0)
                            val *= pcoefficient[ip];
        };

        chain_rule(sx, nxtrue_rdata, nx_rdata, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
        chain_rule(sx0, nxtrue_rdata, nx_rdata, 1);
        chain_rule(sx_ss, nxtrue_rdata, nx_rdata, 1);
    }

    if (o2mode == SecondOrderMode::full) {
        if (!s2llh.empty() && !sllh.empty()) {
            for (int ip = 0; ip < nplist; ++ip) {
                for (int iJ = 1; iJ < nJ; ++iJ) {
                    auto&& val = s2llh[ip * nplist + (iJ - 1)];
                    if (val != 0.0)
                        val *= pcoefficient.at(ip) * augcoefficient[iJ - 1];
                    if (model.plist(ip) == iJ - 1)
                        val += sllh.at(ip) * coefficient.at(ip);
                }
            }
        }

        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty())
                return;

            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));
            Expects(coefficient.size() == gsl::narrow<index_type>(nplist));

            for (int ip = 0; ip < nplist; ++ip)
                for (int iJ = 1; iJ < nJ; ++iJ)
                    for (int i1 = 0; i1 < n1; ++i1)
                        for (int i2 = 0; i2 < n2; ++i2) {
                            auto idx
                                = (i2 * nplist + ip) * stride1 + i1 + iJ * n1;
                            auto&& val = sens.at(idx);
                            if (val != 0.0)
                                val *= pcoefficient.at(ip)
                                       * augcoefficient[iJ - 1];
                            if (model.plist(ip) == iJ - 1)
                                val += sens.at(
                                           (i2 * nplist + ip) * stride1 + i1
                                       )
                                       * coefficient[ip];
                        }
        };

        chain_rule(sx, nxtrue_rdata, nx_rdata, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
    } else if (o2mode == SecondOrderMode::directional) {
        for (int ip = 0; ip < nplist; ++ip) {
            auto&& val = sllh.at(ip);
            if (val != 0.0)
                val *= pcoefficient.at(ip);
            val += model.k()[nk - nplist + ip] * sllh.at(ip)
                   / unscaledParameters[model.plist(ip)];
        }

        auto chain_rule = [&](auto& sens, int n1, int stride1, int n2) {
            if (sens.empty())
                return;

            using index_type =
                typename std::remove_reference_t<decltype(sens)>::size_type;
            Expects(
                sens.size() == gsl::narrow<index_type>(n2 * nplist * stride1)
            );
            Expects(n1 <= stride1);
            Expects(pcoefficient.size() == gsl::narrow<index_type>(nplist));

            for (int ip = 0; ip < nplist; ++ip)
                for (int i1 = 0; i1 < n1; ++i1)
                    for (int i2 = 0; i2 < n2; ++i2) {
                        auto idx = (i2 * nplist + ip) * stride1 + i1 + n1;
                        auto&& val = sens.at(idx);
                        if (val != 0.0)
                            val *= pcoefficient.at(ip);
                        val += model.k()[nk - nplist + ip]
                               * sens.at((i2 * nplist + ip) * stride1 + i1)
                               / unscaledParameters[model.plist(ip)];
                    }
        };

        chain_rule(sx, nxtrue_rdata, nx_rdata, nt);
        chain_rule(sy, nytrue, ny, nt);
        chain_rule(ssigmay, nytrue, ny, nt);
        chain_rule(sz, nztrue, nz, nmaxevent);
        chain_rule(ssigmaz, nztrue, nz, nmaxevent);
        chain_rule(srz, nztrue, nz, nmaxevent);
    }
}

void ReturnData::initialize_objective_function(bool enable_chi2) {
    if (rdata_reporting == RDataReporting::likelihood
        || rdata_reporting == RDataReporting::observables_likelihood
        || rdata_reporting == RDataReporting::full) {
        llh = 0.0;
        std::ranges::fill(sllh, 0.0);
        std::ranges::fill(s2llh, 0.0);
    }
    if ((rdata_reporting == RDataReporting::residuals
         || rdata_reporting == RDataReporting::full)
        && enable_chi2)
        chi2 = 0.0;
}

static realtype
fres(realtype y, realtype my, realtype sigma_y, ObservableScaling scale) {
    switch (scale) {
    case ObservableScaling::lin:
        return (my - y) / sigma_y;
    case ObservableScaling::log:
        return (std::log(my) - std::log(y)) / sigma_y;
    case ObservableScaling::log10:
        return (std::log10(my) - std::log10(y)) / sigma_y;
    default:
        throw std::invalid_argument("only lin, log, log10 allowed.");
    }
}

static realtype fres_error(realtype sigma_y, realtype sigma_offset) {
    return sqrt(2 * log(sigma_y) + sigma_offset);
}

void ReturnData::fres(
    int const it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (res.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.get_observable(y_it, ts[it], sol.x);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.get_observable_sigma(sigmay_it, it, &edata);

    auto observedData = edata.get_measurements_ptr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt = iy + it * edata.nytrue();
        if (!edata.is_set_measurement(it, iy))
            continue;

        res.at(iyt) = amici::fres(
            y_it.at(iy), observedData[iy], sigmay_it.at(iy),
            model.get_observable_scaling(iy)
        );

        if (sigma_res)
            res.at(iyt + nt * nytrue)
                = fres_error(sigmay_it.at(iy), sigma_offset);
    }
}

void ReturnData::fchi2(int const it, ExpData const& edata) {
    if (res.empty() || isnan(chi2))
        return;

    for (int iy = 0; iy < nytrue; ++iy) {
        int iyt_true = iy + it * nytrue;
        chi2 += pow(res.at(iyt_true), 2);
        if (sigma_res && edata.is_set_measurement(it, iy))
            chi2 += pow(res.at(iyt_true + nt * nytrue), 2) - sigma_offset;
    }
}

static realtype fsres(
    realtype y, realtype sy, realtype my, realtype sigma_y, realtype ssigma_y,
    ObservableScaling scale
) {
    auto res = fres(y, my, sigma_y, scale);
    switch (scale) {
    case ObservableScaling::lin:
        return (-sy - ssigma_y * res) / sigma_y;
    case ObservableScaling::log:
        return (-sy / y - ssigma_y * res) / sigma_y;
    case ObservableScaling::log10:
        return (-sy / (y * std::log(10)) - ssigma_y * res) / sigma_y;
    default:
        throw std::invalid_argument("only lin, log, log10 allowed.");
    }
}
static realtype
fsres_error(realtype sigma_y, realtype ssigma_y, realtype sigma_offset) {
    return ssigma_y / (fres_error(sigma_y, sigma_offset) * sigma_y);
}

void ReturnData::fsres(
    int const it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (sres.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.get_observable(y_it, ts[it], sol.x);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.get_observable_sensitivity(sy_it, ts[it], sol.x, sol.sx);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.get_observable_sigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.get_observable_sigma_sensitivity(ssigmay_it, sy_it, it, &edata);

    auto observedData = edata.get_measurements_ptr(it);
    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.is_set_measurement(it, iy))
            continue;
        for (int ip = 0; ip < nplist; ++ip) {
            int idx = (iy + it * edata.nytrue()) * nplist + ip;

            sres.at(idx) = amici::fsres(
                y_it.at(iy), sy_it.at(iy + ny * ip), observedData[iy],
                sigmay_it.at(iy), ssigmay_it.at(iy + ny * ip),
                model.get_observable_scaling(iy)
            );

            if (sigma_res) {
                int idx_res
                    = (iy + it * edata.nytrue() + edata.nytrue() * edata.nt())
                          * nplist
                      + ip;
                sres.at(idx_res) = fsres_error(
                    sigmay_it.at(iy), ssigmay_it.at(iy + ny * ip), sigma_offset
                );
            }
        }
    }
}

void ReturnData::fFIM(
    int it, Model& model, SolutionState const& sol, ExpData const& edata
) {
    if (FIM.empty())
        return;

    std::vector<realtype> y_it(ny, 0.0);
    model.get_observable(y_it, ts[it], sol.x);
    std::vector<realtype> sy_it(ny * nplist, 0.0);
    model.get_observable_sensitivity(sy_it, ts[it], sol.x, sol.sx);

    std::vector<realtype> sigmay_it(ny, 0.0);
    model.get_observable_sigma(sigmay_it, it, &edata);
    std::vector<realtype> ssigmay_it(ny * nplist, 0.0);
    model.get_observable_sigma_sensitivity(ssigmay_it, sy_it, it, &edata);

    /*
     * https://www.wolframalpha.com/input/?i=d%2Fdu+d%2Fdv+log%28s%28u%2Cv%29%29+%2B+0.5+*+%28r%28u%2Cv%29%2Fs%28u%2Cv%29%29%5E2
     * r = (y - m)
     * r_du = y_du
     * d/du(d/(dv)(log(s) + 0.5 (r/s)^2)) =
     * -(2*y_du*s_dv)*r/s^3 - (2*y_dv*s_du)*r/s^3 + y_du_dv*r/s^2
     * 22222222222222222222   2222222222222222222   #############
     *
     * + (y_dv*y_du)/s^2 + (3*s_dv*s_du)*r^2/s^4 - s_du_dv*r^2/s^3
     *   111111111111111   333333333333333333333   ***************
     *
     * - (s_dv*s_du)/s^2 + (s_du_dv)/s
     *                     +++++++++++
     *
     *
     * we compute this using fsres:
     * sres_u * sres_v = (y_du*s - s_du * r) * (y_dv*s - s_dv * r) / s^4
     * = y_du*y_dv/s^2 - (y_du*s_dv + y_dv*s_du)*r/s^3 + s_du*s_dv*r^2/s^4
     *   1111111111111   22222222222222222222222222222   33333333333333333
     *
     * r should be on the same order as s. We keep 1/s^2 and drop 1/s terms.
     * drop:
     * ********: r^2/s^3 term, typically zero anyways
     * ########: r/s^2 term, requires second order sensitivities
     * ++++++++: 1/s term, typically zero anyways
     *
     * keep:
     * 123123123: accounted for by sres*sres,
     * but
     * - (s_dv*s_du)/s^2 is unaccounted for and
     * - (s_dv*y_du + s_du*y_dv)*r/s^3 is missing from 2222 and
     * + 2*(s_du*s_dv)*r^2/s^4 is missing from 3333
     *
     * s_dv*y_du and s_du*y_dv are usually zero since we do not have parameters
     * that affect both observables and sigmas. Accordingly, it is hard to know
     * emprically whether these terms are important or not.
     *
     * This leaves us with
     * + (s_du*s_dv)(2*r^2-s^2)/s^4
     * which may be problematic, since this expression does not factorise and
     * may thus introduce directions of negative curvature.
     *
     * For the least squares trick, where error residuals
     * er = sqrt(log(s) + c), with sensitivity er_du = s_du/(2*s*er). This
     * would yield terms (s_du*s_dv)*(s^2/(4*er^2))/s^4.
     * These terms are guaranteed to yield positive curvature, but go to zero
     * in the limit c -> Infty.
     *
     * Empirically, simply taking this limit and dropping all missing terms,
     * works substantially better. This was evaluated using the fides optimizer
     * on the Boehm2014 Benchmark example.
     */

    auto observedData = edata.get_measurements_ptr(it);

    for (int iy = 0; iy < nytrue; ++iy) {
        if (!edata.is_set_measurement(it, iy))
            continue;
        auto y = y_it.at(iy);
        auto m = observedData[iy];
        auto s = sigmay_it.at(iy);
        auto os = model.get_observable_scaling(iy);
        // auto r = amici::fres(y, m, s);
        for (int ip = 0; ip < nplist; ++ip) {
            auto dy_i = sy_it.at(iy + ny * ip);
            auto ds_i = ssigmay_it.at(iy + ny * ip);
            auto sr_i = amici::fsres(y, dy_i, m, s, ds_i, os);
            realtype sre_i = 0.0;
            if (sigma_res)
                sre_i = fsres_error(s, ds_i, sigma_offset);
            for (int jp = 0; jp < nplist; ++jp) {
                auto dy_j = sy_it.at(iy + ny * jp);
                auto ds_j = ssigmay_it.at(iy + ny * jp);
                auto sr_j = amici::fsres(y, dy_j, m, s, ds_j, os);
                FIM.at(ip + nplist * jp) += sr_i * sr_j;
                if (sigma_res) {
                    auto sre_j = fsres_error(s, ds_j, sigma_offset);
                    FIM.at(ip + nplist * jp) += sre_i * sre_j;
                }
                /*+ ds_i*ds_j*(2*pow(r/pow(s,2.0), 2.0) - 1/pow(s,2.0));*/
            }
        }
    }
}

ModelContext::ModelContext(Model* model)
    : model_(model)
    , original_state_(model->get_model_state()) {}

ModelContext::~ModelContext() noexcept(false) { restore(); }

void ModelContext::restore() { model_->set_model_state(original_state_); }

} // namespace amici
