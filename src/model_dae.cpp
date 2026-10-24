#include "amici/model_dae.h"
#include "amici/solver_idas.h"
#include "amici/sundials_linsol_wrapper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace amici {

void Model_DAE::fJ(
    realtype const t, realtype const cj, AmiVector const& x,
    AmiVector const& dx, AmiVector const& xdot, SUNMatrix J
) {
    fJ(t, cj, x.get_nvector(), dx.get_nvector(), xdot.get_nvector(), J);
}

void Model_DAE::fJ(
    realtype const t, realtype const cj, const_N_Vector x, const_N_Vector dx,
    const_N_Vector /*xdot*/, SUNMatrix J
) {
    fJSparse(t, cj, x, dx, derived_state_.J_);
    derived_state_.J_.refresh();
    auto JDense = SUNMatrixWrapper(J);
    derived_state_.J_.to_dense(JDense);
}

void Model_DAE::fJSparse(
    realtype const t, realtype const cj, AmiVector const& x,
    AmiVector const& dx, AmiVector const& /*xdot*/, SUNMatrix J
) {
    fJSparse(t, cj, x.get_nvector(), dx.get_nvector(), J);
}

void Model_DAE::fJSparse(
    realtype const t, realtype const cj, const_N_Vector x, const_N_Vector dx,
    SUNMatrix J
) {
    auto const x_pos = compute_x_pos(x);
    fdwdx(t, N_VGetArrayPointerConst(x_pos), false);

    auto JSparse = SUNMatrixWrapper(J);
    // python generated
    derived_state_.dxdotdx_explicit.zero();
    derived_state_.dxdotdx_implicit.zero();
    if (derived_state_.dxdotdx_explicit.capacity()) {
        fdxdotdx_explicit_colptrs(derived_state_.dxdotdx_explicit);
        fdxdotdx_explicit_rowvals(derived_state_.dxdotdx_explicit);
        fdxdotdx_explicit(
            derived_state_.dxdotdx_explicit.data(), t,
            N_VGetArrayPointerConst(x_pos), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), state_.h.data(),
            N_VGetArrayPointerConst(dx), derived_state_.w_.data()
        );
    }
    fdxdotdw(t, x_pos, dx);
    /* Sparse matrix multiplication
     dxdotdx_implicit += dxdotdw * dwdx */
    derived_state_.dxdotdw_.sparse_multiply(
        derived_state_.dxdotdx_implicit, derived_state_.dwdx_
    );

    derived_state_.dfdx_.sparse_add(
        derived_state_.dxdotdx_explicit, 1.0, derived_state_.dxdotdx_implicit,
        1.0
    );
    fM(t, x_pos);
    JSparse.sparse_add(derived_state_.MSparse_, -cj, derived_state_.dfdx_, 1.0);
}

void Model_DAE::fJv(
    realtype const t, AmiVector const& x, AmiVector const& dx,
    AmiVector const& /*xdot*/, AmiVector const& v, AmiVector& Jv,
    realtype const cj
) {
    fJv(t, x.get_nvector(), dx.get_nvector(), v.get_nvector(), Jv.get_nvector(),
        cj);
}

void Model_DAE::fJv(
    realtype const t, const_N_Vector x, const_N_Vector dx, const_N_Vector v,
    N_Vector Jv, realtype const cj
) {
    N_VConst(0.0, Jv);
    fJSparse(t, cj, x, dx, derived_state_.J_);
    derived_state_.J_.refresh();
    derived_state_.J_.multiply(Jv, v);
}

void Model_DAE::froot(
    realtype const t, AmiVector const& x, AmiVector const& dx,
    gsl::span<realtype> const root
) {
    froot(t, x.get_nvector(), dx.get_nvector(), root);
}

void Model_DAE::froot(
    realtype const t, const_N_Vector x, const_N_Vector dx,
    gsl::span<realtype> root
) {
    std::ranges::fill(root, 0.0);
    auto const x_pos = compute_x_pos(x);
    // TODO(performance) only dynamic expressions? consider storing a flag
    //  for whether froot actually depends on `w`.
    fw(t, N_VGetArrayPointerConst(x_pos));
    froot(
        root.data(), t, N_VGetArrayPointerConst(x_pos),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data(), derived_state_.w_.data(), N_VGetArrayPointerConst(dx)
    );
}

void Model_DAE::fxdot(
    realtype const t, AmiVector const& x, AmiVector const& dx, AmiVector& xdot
) {
    fxdot(t, x.get_nvector(), dx.get_nvector(), xdot.get_nvector());
}

void Model_DAE::fxdot(
    realtype const t, const_N_Vector x, const_N_Vector dx, N_Vector xdot
) {
    auto const x_pos = compute_x_pos(x);
    fw(t, N_VGetArrayPointerConst(x), false);
    N_VConst(0.0, xdot);
    fxdot(
        N_VGetArrayPointer(xdot), t, N_VGetArrayPointerConst(x_pos),
        state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
        state_.h.data(), N_VGetArrayPointerConst(dx), derived_state_.w_.data()
    );
}

void Model_DAE::fdxdt(
    realtype const t, AmiVector const& x, AmiVector const& dx, AmiVector& dxdt
) {
    // residual of the differential equations: f(x) - dx
    fxdot(t, x, dx, dxdt);
    for (int ix = 0; ix < nx_solver; ++ix) {
        if (get_id_list().at(ix) == 1.0)
            dxdt[ix] += dx[ix];
        else
            dxdt[ix] = dx[ix];
    }
}

void Model_DAE::fJDiag(
    realtype const t, AmiVector& JDiag, realtype const /*cj*/,
    AmiVector const& x, AmiVector const& dx
) {
    fJSparse(t, 0.0, x.get_nvector(), dx.get_nvector(), derived_state_.J_);
    derived_state_.J_.refresh();
    derived_state_.J_.to_diag(JDiag.get_nvector());
    if (check_finite(JDiag.get_vector(), ModelQuantity::JDiag, t)
        != AMICI_SUCCESS)
        throw AmiException("Evaluation of fJDiag failed!");
}

void Model_DAE::fdxdotdw(
    realtype const t, const_N_Vector x, const_N_Vector const dx
) {
    derived_state_.dxdotdw_.zero();
    if (nw > 0 && derived_state_.dxdotdw_.capacity()) {
        auto const x_pos = compute_x_pos(x);

        fdxdotdw_colptrs(derived_state_.dxdotdw_);
        fdxdotdw_rowvals(derived_state_.dxdotdw_);
        fdxdotdw(
            derived_state_.dxdotdw_.data(), t, N_VGetArrayPointerConst(x_pos),
            state_.unscaled_parameters.data(), state_.fixed_parameters.data(),
            state_.h.data(), N_VGetArrayPointerConst(dx),
            derived_state_.w_.data()
        );
    }
}

void Model_DAE::fdxdotdp(
    realtype const t, const_N_Vector const x, const_N_Vector const dx
) {
    auto const x_pos = compute_x_pos(x);
    fdwdp(t, N_VGetArrayPointerConst(x_pos));

    derived_state_.dxdotdp_explicit.zero();
    derived_state_.dxdotdp_implicit.zero();
    if (derived_state_.dxdotdp_explicit.capacity()) {
        fdxdotdp_explicit_colptrs(derived_state_.dxdotdp_explicit);
        fdxdotdp_explicit_rowvals(derived_state_.dxdotdp_explicit);
        fdxdotdp_explicit(
            derived_state_.dxdotdp_explicit.data(), t,
            N_VGetArrayPointerConst(x_pos), state_.unscaled_parameters.data(),
            state_.fixed_parameters.data(), state_.h.data(),
            N_VGetArrayPointerConst(dx), derived_state_.w_.data()
        );
    }

    fdxdotdw(t, x_pos, dx);
    /* Sparse matrix multiplication
     dxdotdp_implicit += dxdotdw * dwdp */
    derived_state_.dxdotdw_.sparse_multiply(
        derived_state_.dxdotdp_implicit, derived_state_.dwdp_
    );

    derived_state_.dxdotdp_full.sparse_add(
        derived_state_.dxdotdp_explicit, 1.0, derived_state_.dxdotdp_implicit,
        1.0
    );
}

void Model_DAE::initialize_state_derivatives(
    realtype const t, AmiVector& x, AmiVector& dx, AmiVectorArray const& sx,
    AmiVectorArray& sdx, bool const computeSensitivities
) {
    auto const& id = get_id_list();
    // For the differential states (mass matrix entry 1), the residual for
    // zero derivatives is the right-hand side.
    dx.zero();
    AmiVector residual(dx);
    fxdot(t, x, dx, residual);
    for (int ix = 0; ix < nx_solver; ++ix)
        dx[ix] = id[ix] != 0.0 ? residual[ix] : 0.0;

    if (!computeSensitivities)
        return;

    AmiVector const zero(nx_solver, dx.get_ctx());
    for (int ip = 0; ip < nplist(); ++ip) {
        fsxdot(t, x, dx, ip, sx[ip], zero, residual);
        for (int ix = 0; ix < nx_solver; ++ix)
            sdx[ip][ix] = id[ix] != 0.0 ? residual[ix] : 0.0;
    }
}

namespace {
/**
 * @brief Indices of the algebraic states, and the position of each state
 * among them (-1 for differential states).
 */
std::pair<std::vector<int>, std::vector<int>>
algebraic_state_indices(std::vector<realtype> const& id) {
    std::vector<int> alg_states;
    std::vector<int> alg_pos(id.size(), -1);
    for (std::size_t ix = 0; ix < id.size(); ++ix) {
        if (id[ix] == 0.0) {
            alg_pos[ix] = gsl::narrow<int>(alg_states.size());
            alg_states.push_back(gsl::narrow<int>(ix));
        }
    }
    return {alg_states, alg_pos};
}
} // namespace

void Model_DAE::setup_algebraic_jacobian_solver(
    realtype const t, AmiVector const& x, std::vector<int> const& alg_pos,
    bool const transpose, char const* purpose, SUNLinSolDense& linsol
) {
    // the algebraic equations do not depend on dx
    AmiVector const zero(nx_solver, derived_state_.sunctx_);
    fJSparse(t, 0.0, x.get_nvector(), zero.get_nvector(), derived_state_.J_);
    derived_state_.J_.refresh();
    auto const& J = derived_state_.J_;

    auto& A = linsol.getMatrix();
    A.zero();
    for (int ix = 0; ix < nx_solver; ++ix) {
        auto const col = alg_pos[ix];
        if (col < 0)
            continue;
        for (auto idx = J.get_indexptr(ix); idx < J.get_indexptr(ix + 1);
             ++idx) {
            if (auto const row = alg_pos[J.get_indexval(idx)]; row >= 0) {
                if (transpose)
                    A.set_data(col, row, J.get_data(idx));
                else
                    A.set_data(row, col, J.get_data(idx));
            }
        }
    }
    try {
        linsol.setup();
    } catch (AmiException const&) {
        throw AmiException(
            "Singular algebraic Jacobian at t = %g: cannot %s.", t, purpose
        );
    }
}

void Model_DAE::project_algebraic_objective_derivatives(
    realtype const t, AmiVector const& x, gsl::span<realtype> dJydx,
    gsl::span<realtype> dJydp
) {
    // With differential states x_d and algebraic states x_a determined by
    // 0 = f_a(x_d, x_a, p), the objective J(x_d, x_a, p) has the total
    // derivatives
    //   dJ/dx_d = J_{x_d} - (df_a/dx_d)^T mu,
    //   dJ/dp   = J_p     - (df_a/dp)^T mu,
    // where (df_a/dx_a)^T mu = J_{x_a}.
    // Algebraic equation i determines algebraic state i. Only first-order
    // derivatives (entries iJ == 0, stride nJ) are handled.
    auto const [alg_states, alg_pos] = algebraic_state_indices(get_id_list());
    auto const n_alg = gsl::narrow<int>(alg_states.size());
    if (n_alg == 0)
        return;

    AmiVector mu(n_alg, derived_state_.sunctx_);
    bool any_nonzero = false;
    for (int ia = 0; ia < n_alg; ++ia) {
        mu[ia] = dJydx[alg_states[ia] * nJ];
        any_nonzero = any_nonzero || mu[ia] != 0.0;
    }
    if (!any_nonzero)
        return;

    // solve (df_a/dx_a)^T mu = J_{x_a}
    SUNLinSolDense linsol(mu);
    setup_algebraic_jacobian_solver(
        t, x, alg_pos, true,
        "map objective derivatives w.r.t. algebraic states", linsol
    );
    auto const& J = derived_state_.J_;
    if (auto const status
        = linsol.solve(mu.get_nvector(), mu.get_nvector(), 0.0);
        status != SUN_SUCCESS) {
        throw AmiException(
            "Mapping objective derivatives w.r.t. algebraic states failed "
            "with code %d",
            status
        );
    }

    // mu as a full-length vector (zero for the differential equations)
    std::vector<realtype> mu_full(nx_solver, 0.0);
    for (int ia = 0; ia < n_alg; ++ia)
        mu_full[alg_states[ia]] = mu[ia];

    // dJ/dx_d -= (df_a/dx_d)^T mu; dJ/dx_a is now accounted for via mu
    std::vector<int> all_states(nx_solver);
    std::iota(all_states.begin(), all_states.end(), 0);
    std::vector<realtype> correction(nx_solver, 0.0);
    J.multiply(correction, mu_full, all_states, true);
    for (int ix = 0; ix < nx_solver; ++ix) {
        dJydx[ix * nJ]
            = alg_pos[ix] >= 0 ? 0.0 : dJydx[ix * nJ] - correction[ix];
    }

    if (dJydp.empty())
        return;

    // dJ/dp -= (df_a/dp)^T mu; the algebraic equations do not depend on dx
    AmiVector const zero(nx_solver, derived_state_.sunctx_);
    fdxdotdp(t, x.get_nvector(), zero.get_nvector());
    correction.assign(nplist(), 0.0);
    derived_state_.dxdotdp_full.multiply(
        correction, mu_full, state_.plist, true
    );
    for (int ip = 0; ip < nplist(); ++ip)
        dJydp[ip * nJ] -= correction[ip];
}

void Model_DAE::make_algebraic_derivatives_consistent(
    realtype const t, AmiVector const& x, AmiVector& dx, AmiVectorArray* sx
) {
    // For the algebraic equations 0 = f_a(t, x, p), the time derivatives and
    // the sensitivities v of the states satisfy
    //   0 = df_a/dx v + r,
    // with r = df_a/dt or r = df_a/dp, respectively. The algebraic
    // components of v are corrected such that this holds, which is a single
    // Newton step since the equations are linear in v.
    // Algebraic equation i determines algebraic state i.
    auto const [alg_states, alg_pos] = algebraic_state_indices(get_id_list());
    auto const n_alg = gsl::narrow<int>(alg_states.size());
    if (n_alg == 0)
        return;

    AmiVector const zero(nx_solver, derived_state_.sunctx_);
    AmiVector residual(nx_solver, derived_state_.sunctx_);
    fdxdotdt(t, x, residual);

    AmiVector b(n_alg, derived_state_.sunctx_);
    SUNLinSolDense linsol(b);
    setup_algebraic_jacobian_solver(
        t, x, alg_pos, false,
        "compute consistent derivatives of the algebraic states", linsol
    );
    auto const& J = derived_state_.J_;

    // v_a -= (df_a/dx_a)^{-1} residual_a
    auto correct = [&](AmiVector& v, AmiVector const& r) {
        for (int ia = 0; ia < n_alg; ++ia)
            b[ia] = r[alg_states[ia]];
        if (auto const status
            = linsol.solve(b.get_nvector(), b.get_nvector(), 0.0);
            status != SUN_SUCCESS) {
            throw AmiException(
                "Computing consistent derivatives of the algebraic states "
                "failed with code %d",
                status
            );
        }
        for (int ia = 0; ia < n_alg; ++ia)
            v[alg_states[ia]] -= b[ia];
    };

    // `residual` holds df_a/dt
    J.multiply(residual.get_nvector(), dx.get_nvector());
    correct(dx, residual);

    if (!sx)
        return;
    for (int ip = 0; ip < nplist(); ++ip) {
        // df_a/dp + df_a/dx sx
        fsxdot(t, x, zero, ip, (*sx)[ip], zero, residual);
        correct((*sx)[ip], residual);
    }
}

void Model_DAE::fdxdotdt(
    realtype const t, AmiVector const& x, AmiVector& dxdotdt
) {
    auto const x_pos = compute_x_pos(x.get_nvector());
    auto const* x_data = N_VGetArrayPointerConst(x_pos);
    AmiVector const zero(nx_solver, derived_state_.sunctx_);
    // w and dwdw
    fdwdx(t, x_data, false);

    dxdotdt.zero();
    fdxdotdt_explicit(
        dxdotdt.data(), t, x_data, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data(), state_.h.data(), zero.data(),
        derived_state_.w_.data()
    );
    if (!nw)
        return;

    // dxdotdt += dxdotdw * dwdt, with the total time derivative of w at
    // fixed x: dwdt = sum_k dwdw^k dwdt_explicit
    // Time derivatives of splines are not implemented yet (they would have to
    // be evaluated here). NaN makes any algebraic equation that depends on
    // them non-finite, which is caught below.
    std::vector<realtype> const dspl(
        nspl, std::numeric_limits<realtype>::quiet_NaN()
    );
    std::vector<realtype> dwdt_k(nw, 0.0);
    fdwdt(
        dwdt_k.data(), t, x_data, state_.unscaled_parameters.data(),
        state_.fixed_parameters.data(), state_.h.data(),
        derived_state_.w_.data(), state_.total_cl.data(),
        derived_state_.spl_.data(), dspl.data()
    );
    auto dwdt = dwdt_k;
    std::vector<realtype> tmp(nw);
    for (int irecursion = 1; irecursion <= w_recursion_depth; irecursion++) {
        std::ranges::fill(tmp, 0.0);
        derived_state_.dwdw_.multiply(tmp, dwdt_k);
        std::swap(tmp, dwdt_k);
        for (int iw = 0; iw < nw; ++iw)
            dwdt[iw] += dwdt_k[iw];
    }
    fdxdotdw(t, x_pos, zero.get_nvector());
    derived_state_.dxdotdw_.multiply(
        gsl::make_span(dxdotdt.data(), nx_solver), dwdt
    );
    auto const& id = get_id_list();
    for (int ix = 0; ix < nx_solver; ++ix) {
        if (id[ix] != 0.0)
            dxdotdt[ix] = 0.0;
        else if (!std::isfinite(dxdotdt[ix]))
            throw AmiException(
                "Non-finite time derivative of algebraic equation %d at "
                "t = %g. Note that time derivatives of splines are not "
                "supported yet.",
                ix, t
            );
    }
}

void Model_DAE::fM(realtype const /*t*/, const_N_Vector /*x*/) {
    /*
     * non-algebraic states in python generated code always have factor
     * 1 in the mass matrix, so we can easily construct this matrix here
     * and avoid having to generate c++ code
     */
    int ndiff = 0;
    for (int ix = 0; ix < nx_solver; ix++) {
        derived_state_.MSparse_.set_indexptr(ix, ndiff);
        if (get_id_list().at(ix) == 1.0) {
            derived_state_.MSparse_.set_data(ndiff, 1.0);
            derived_state_.MSparse_.set_indexval(ndiff, ix);
            ndiff++;
        }
    }
    derived_state_.MSparse_.set_indexptr(nx_solver, ndiff);
    assert(ndiff == derived_state_.MSparse_.capacity());
}

std::unique_ptr<Solver> Model_DAE::create_solver() {
    return std::unique_ptr<Solver>(new IDASolver());
}

void Model_DAE::fJSparse(
    SUNMatrixContent_Sparse /*JSparse*/, realtype /*t*/, realtype const* /*x*/,
    double const* /*p*/, double const* /*k*/, realtype const* /*h*/,
    realtype /*cj*/, realtype const* /*dx*/, realtype const* /*w*/,
    realtype const* /*dwdx*/
) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::froot(
    realtype* /*root*/, realtype const /*t*/, realtype const* /*x*/,
    double const* /*p*/, double const* /*k*/, realtype const* /*h*/,
    realtype const* /*w*/, realtype const* /*dx*/
) {
    throw AmiException(
        "Requested functionality is not supported as %s is not "
        "implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdp_explicit(
    realtype* /*dxdotdp_explicit*/, realtype const /*t*/, realtype const* /*x*/,
    realtype const* /*p*/, realtype const* /*k*/, realtype const* /*h*/,
    realtype const* /*dx*/, realtype const* /*w*/
) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdp_explicit_colptrs(SUNMatrixWrapper& /*dxdotdp*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdp_explicit_rowvals(SUNMatrixWrapper& /*dxdotdp*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdx_explicit(
    realtype* /*dxdotdx_explicit*/, realtype const /*t*/, realtype const* /*x*/,
    realtype const* /*p*/, realtype const* /*k*/, realtype const* /*h*/,
    realtype const* /*dx*/, realtype const* /*w*/
) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdx_explicit_colptrs(SUNMatrixWrapper& /*dxdotdx*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdx_explicit_rowvals(SUNMatrixWrapper& /*dxdotdx*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdw(
    realtype* /*dxdotdw*/, realtype const /*t*/, realtype const* /*x*/,
    realtype const* /*p*/, realtype const* /*k*/, realtype const* /*h*/,
    realtype const* /*dx*/, realtype const* /*w*/
) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdw_colptrs(SUNMatrixWrapper& /*dxdotdw*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::fdxdotdw_rowvals(SUNMatrixWrapper& /*dxdotdw*/) {
    throw AmiException(
        "Requested functionality is not supported as %s "
        "is not implemented for this model!",
        __func__
    ); // not implemented
}

void Model_DAE::
    fM(realtype* /*M*/, realtype const /*t*/, realtype const* /*x*/,
       realtype const* /*p*/, realtype const* /*k*/) {}

void Model_DAE::fJB(
    realtype const t, realtype cj, AmiVector const& x, AmiVector const& dx,
    AmiVector const& xB, AmiVector const& /*dxB*/, AmiVector const& /*xBdot*/,
    SUNMatrix JB
) {
    fJB(t, cj, x.get_nvector(), dx.get_nvector(), xB.get_nvector(),
        dx.get_nvector(), JB);
}

void Model_DAE::fJB(
    realtype const t, realtype const cj, const_N_Vector x, const_N_Vector dx,
    const_N_Vector /*xB*/, const_N_Vector /*dxB*/, SUNMatrix JB
) {
    // JB = -(df/dx)^T - cj * M^T = -(J(-cj))^T, see fxBdot
    fJSparse(t, -cj, x, dx, derived_state_.J_);
    derived_state_.J_.refresh();
    auto JBDense = SUNMatrixWrapper(JB);
    derived_state_.J_.transpose(JBDense, -1.0, nxtrue_solver);
}

void Model_DAE::fJSparseB(
    realtype const t, realtype const cj, AmiVector const& x,
    AmiVector const& dx, AmiVector const& xB, AmiVector const& dxB,
    AmiVector const& /*xBdot*/, SUNMatrix JB
) {
    fJSparseB(
        t, cj, x.get_nvector(), dx.get_nvector(), xB.get_nvector(),
        dxB.get_nvector(), JB
    );
}

void Model_DAE::fJSparseB(
    realtype const t, realtype const cj, const_N_Vector x, const_N_Vector dx,
    const_N_Vector /*xB*/, const_N_Vector /*dxB*/, SUNMatrix JB
) {
    // JB = -(df/dx)^T - cj * M^T = -(J(-cj))^T, see fxBdot
    fJSparse(t, -cj, x, dx, derived_state_.J_);
    derived_state_.J_.refresh();
    auto JSparseB = SUNMatrixWrapper(JB);
    derived_state_.J_.transpose(JSparseB, -1.0, nxtrue_solver);
}

void Model_DAE::fJvB(
    realtype const t, const_N_Vector x, const_N_Vector dx, const_N_Vector xB,
    const_N_Vector dxB, const_N_Vector vB, N_Vector JvB, realtype cj
) {
    N_VConst(0.0, JvB);
    fJSparseB(t, cj, x, dx, xB, dxB, derived_state_.JB_);
    derived_state_.JB_.refresh();
    derived_state_.JB_.multiply(JvB, vB);
}

void Model_DAE::fxBdot(
    realtype const t, const_N_Vector x, const_N_Vector dx, const_N_Vector xB,
    const_N_Vector dxB, N_Vector xBdot
) {
    // Residual of the adjoint DAE for F = f(x) - M * dx:
    //   0 = -M^T * dxB - (df/dx)^T * xB
    N_VConst(0.0, xBdot);
    fJSparseB(t, 0.0, x, dx, xB, dxB, derived_state_.JB_);
    derived_state_.JB_.refresh();
    derived_state_.JB_.multiply(xBdot, xB);
    fM(t, x);
    derived_state_.MSparse_.multiply(xBdot, dxB, -1.0);
}

void Model_DAE::fqBdot(
    realtype const t, const_N_Vector x, const_N_Vector dx, const_N_Vector xB,
    const_N_Vector /*dxB*/, N_Vector qBdot
) {
    N_VConst(0.0, qBdot);
    fdxdotdp(t, x, dx);
    derived_state_.dxdotdp_full.multiply(qBdot, xB, state_.plist, true);
    N_VScale(-1.0, qBdot, qBdot);
}

void Model_DAE::fxBdot_ss(
    realtype const t, AmiVector const& xB, AmiVector const& dxB,
    AmiVector& xBdot
) {
    fxBdot_ss(t, xB.get_nvector(), dxB.get_nvector(), xBdot.get_nvector());
}

void Model_DAE::fxBdot_ss(
    realtype /*t*/, const_N_Vector xB, const_N_Vector /*dxB*/, N_Vector xBdot
) const {
    /* Right hand side of the adjoint state for steady state computations.
     J is fixed (as x remains in steady state), so the RHS becomes simple. */
    N_VConst(0.0, xBdot);
    derived_state_.JB_.multiply(xBdot, xB);
    /* Mind the minus sign... */
    N_VScale(-1.0, xBdot, xBdot);
}

void Model_DAE::fqBdot_ss(
    realtype /*t*/, const_N_Vector xB, const_N_Vector /*dxB*/, N_Vector qBdot
) const {
    /* Quadratures when computing adjoints for steady state. The integrand is
     just the adjoint state itself. */
    N_VScale(1.0, const_cast<N_Vector>(xB), qBdot);
}

void Model_DAE::fJSparseB_ss(SUNMatrix JB) {
    /* Just pass the model Jacobian on to JB */
    SUNMatCopy(derived_state_.JB_, JB);
    derived_state_.JB_.refresh();
}

void Model_DAE::write_steady_state_JB(
    realtype const t, realtype cj, AmiVector const& x, AmiVector const& dx,
    AmiVector const& xB, AmiVector const& dxB, AmiVector const& /*xBdot*/
) {
    /* Get backward Jacobian */
    fJSparseB(
        t, cj, x.get_nvector(), dx.get_nvector(), xB.get_nvector(),
        dxB.get_nvector(), derived_state_.JB_
    );
    derived_state_.JB_.refresh();
    /* Switch sign, as we integrate forward in time, not backward */
    derived_state_.JB_.scale(-1);
}

void Model_DAE::fsxdot(
    realtype const t, AmiVector const& x, AmiVector const& dx, int const ip,
    AmiVector const& sx, AmiVector const& sdx, AmiVector& sxdot
) {
    fsxdot(
        t, x.get_nvector(), dx.get_nvector(), ip, sx.get_nvector(),
        sdx.get_nvector(), sxdot.get_nvector()
    );
}

void Model_DAE::fsxdot(
    realtype t, const_N_Vector x, const_N_Vector dx, int ip, const_N_Vector sx,
    const_N_Vector sdx, N_Vector sxdot
) {
    if (ip == 0) {
        // we only need to call this for the first parameter index will be
        // the same for all remaining
        fM(t, x);
        fdxdotdp(t, x, dx);
        fJSparse(t, 0.0, x, dx, derived_state_.J_);
        derived_state_.J_.refresh();
    }

    /* copy dxdotdp and the implicit version over */
    // initialize
    N_VConst(0.0, sxdot);
    realtype* sxdot_tmp = N_VGetArrayPointer(sxdot);
    derived_state_.dxdotdp_full.scatter(
        plist(ip), 1.0, nullptr, gsl::make_span(sxdot_tmp, nx_solver), 0,
        nullptr, 0
    );

    derived_state_.J_.multiply(sxdot, sx);
    derived_state_.MSparse_.multiply(sxdot, sdx, -1.0);
}

} // namespace amici
