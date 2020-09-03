#include <amici/sundials_matrix_wrapper.h>
#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

namespace amici {

void Model_ODE::fJ(const realtype t, const realtype /*cj*/, const AmiVector &x,
                   const AmiVector & /*dx*/, const AmiVector &xdot,
                   SUNMatrix J) {
    fJ(t, x.getNVector(), xdot.getNVector(), J);
}

void Model_ODE::fJ(realtype t, N_Vector x, N_Vector /*xdot*/, SUNMatrix J) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointer(x_pos));
    fJSparse(t, x, J_.get());
    J_.to_dense(J);
}

void Model_ODE::fJSparse(const realtype t, const realtype /*cj*/,
                         const AmiVector &x, const AmiVector & /*dx*/,
                         const AmiVector & /*xdot*/, SUNMatrix J) {
    fJSparse(t, x.getNVector(), J);
}

void Model_ODE::fJSparse(realtype t, N_Vector x, SUNMatrix J) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointer(x_pos));
    SUNMatZero(J);
    if (pythonGenerated) {
        fJSparse(SM_DATA_S(J), t, N_VGetArrayPointer(x_pos),
                 state_.unscaledParameters.data(), state_.fixedParameters.data(),
                 state_.h.data(), w_.data(), dwdx_.data());
        fJSparse_colptrs(SM_INDEXPTRS_S(J));
        fJSparse_rowvals(SM_INDEXVALS_S(J));
    } else {
        fJSparse(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(J)), t,
                 N_VGetArrayPointer(x_pos), state_.unscaledParameters.data(),
                 state_.fixedParameters.data(), state_.h.data(), w_.data(),
                 dwdx_.data());
    }
}

void Model_ODE::fJv(const realtype t, const AmiVector &x,
                    const AmiVector & /*dx*/, const AmiVector & /*xdot*/,
                    const AmiVector &v, AmiVector &Jv, const realtype /*cj*/) {
    fJv(v.getNVector(), Jv.getNVector(), t, x.getNVector());
}

void Model_ODE::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x) {
    N_VConst(0.0, Jv);
    fJSparse(t, x, J_.get());
    J_.multiply(Jv, v);
}

void Model_ODE::froot(const realtype t, const AmiVector &x,
                      const AmiVector & /*dx*/, gsl::span<realtype> root) {
    froot(t, x.getNVector(), root);
}

void Model_ODE::froot(realtype t, N_Vector x, gsl::span<realtype> root) {
    auto x_pos = computeX_pos(x);
    std::fill(root.begin(), root.end(), 0.0);
    froot(root.data(), t, N_VGetArrayPointer(x_pos),
          state_.unscaledParameters.data(), state_.fixedParameters.data(),
          state_.h.data());
}

void Model_ODE::fxdot(const realtype t, const AmiVector &x,
                      const AmiVector & /*dx*/, AmiVector &xdot) {
    fxdot(t, x.getNVector(), xdot.getNVector());
}

void Model_ODE::fxdot(realtype t, N_Vector x, N_Vector xdot) {
    auto x_pos = computeX_pos(x);
    fw(t, N_VGetArrayPointer(x_pos));
    N_VConst(0.0, xdot);
    fxdot(N_VGetArrayPointer(xdot), t, N_VGetArrayPointer(x_pos),
          state_.unscaledParameters.data(), state_.fixedParameters.data(),
          state_.h.data(), w_.data());
}

void Model_ODE::fJDiag(const realtype t, AmiVector &JDiag,
                       const realtype /*cj*/, const AmiVector &x,
                       const AmiVector & /*dx*/) {
    fJDiag(t, JDiag.getNVector(), x.getNVector());
    if (checkFinite(JDiag.getVector(), "Jacobian") != AMICI_SUCCESS)
        throw AmiException("Evaluation of fJDiag failed!");
}

void Model_ODE::fdxdotdw(const realtype t, const N_Vector x) {
    if (nw > 0 && ndxdotdw > 0) {
        auto x_pos = computeX_pos(x);
        dxdotdw_.reset();
        fdxdotdw_colptrs(dxdotdw_.indexptrs());
        fdxdotdw_rowvals(dxdotdw_.indexvals());
        fdxdotdw(dxdotdw_.data(), t, N_VGetArrayPointer(x_pos),
                 state_.unscaledParameters.data(), state_.fixedParameters.data(),
                 state_.h.data(), w_.data());
    }
}

void Model_ODE::fdxdotdp(const realtype t, const N_Vector x) {
    auto x_pos = computeX_pos(x);
    fdwdp(t, N_VGetArrayPointer(x_pos));
    fdxdotdw(t, x_pos);

    if (pythonGenerated) {
        // python generated
        if (dxdotdp_explicit.capacity()) {
            dxdotdp_explicit.reset();
            fdxdotdp_explicit_colptrs(dxdotdp_explicit.indexptrs());
            fdxdotdp_explicit_rowvals(dxdotdp_explicit.indexvals());
            fdxdotdp_explicit(
                dxdotdp_explicit.data(), t, N_VGetArrayPointer(x_pos),
                state_.unscaledParameters.data(), state_.fixedParameters.data(),
                state_.h.data(), w_.data());
        }
        if (nw > 0) {
            /* Sparse matrix multiplication
             dxdotdp_implicit += dxdotdw * dwdp */
            dxdotdp_implicit.reset();
            dxdotdw_.sparse_multiply(dxdotdp_implicit, dwdp_);
        }
        dxdotdp_full.sparse_add(dxdotdp_explicit, 1.0, dxdotdp_implicit, 1.0);
    } else {
        // matlab generated
        for (int ip = 0; ip < nplist(); ip++) {
            N_VConst(0.0, dxdotdp.getNVector(ip));
            fdxdotdp(dxdotdp.data(ip), t, N_VGetArrayPointer(x_pos),
                     state_.unscaledParameters.data(),
                     state_.fixedParameters.data(), state_.h.data(), plist(ip),
                     w_.data(), dwdp_.data());
        }
    }
}

void Model_ODE::fdxdotdp(const realtype t, const AmiVector &x,
                         const AmiVector & /*dx*/) {
    fdxdotdp(t, x.getNVector());
}

std::unique_ptr<Solver> Model_ODE::getSolver() {
    return std::unique_ptr<Solver>(new amici::CVodeSolver());
}

void Model_ODE::fJSparse(SUNMatrixContent_Sparse /*JSparse*/,
                         const realtype /*t*/, const realtype * /*x*/,
                         const realtype * /*p*/, const realtype * /*k*/,
                         const realtype * /*h*/, const realtype * /*w*/,
                         const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparse(realtype * /*JSparse*/, const realtype /*t*/,
                         const realtype * /*x*/, const realtype * /*p*/,
                         const realtype * /*k*/, const realtype * /*h*/,
                         const realtype * /*w*/, const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparse_colptrs(sunindextype * /*indexptrs*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparse_rowvals(sunindextype * /*indexvals*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::froot(realtype * /*root*/, const realtype /*t*/,
                      const realtype * /*x*/, const realtype * /*p*/,
                      const realtype * /*k*/, const realtype * /*h*/) {
    throw AmiException("Requested functionality is not supported as %s is not "
                       "implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdp(realtype * /*dxdotdp*/, const realtype /*t*/,
                         const realtype * /*x*/, const realtype * /*p*/,
                         const realtype * /*k*/, const realtype * /*h*/,
                         const int /*ip*/, const realtype * /*w*/,
                         const realtype * /*dwdp*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdp_explicit(realtype * /*dxdotdp_explicit*/, const realtype /*t*/,
                                  const realtype * /*x*/, const realtype * /*p*/,
                                  const realtype * /*k*/, const realtype * /*h*/,
                                  const realtype * /*w*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdp_explicit_colptrs(sunindextype * /*indexptrs*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdp_explicit_rowvals(sunindextype * /*indexvals*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdw(realtype * /*dxdotdw*/, const realtype /*t*/,
                         const realtype * /*x*/, const realtype * /*p*/,
                         const realtype * /*k*/, const realtype * /*h*/,
                         const realtype * /*w*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdw_colptrs(sunindextype * /*indexptrs*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fdxdotdw_rowvals(sunindextype * /*indexvals*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJB(const realtype t, realtype /*cj*/, const AmiVector &x,
                     const AmiVector & /*dx*/, const AmiVector &xB,
                     const AmiVector & /*dxB*/, const AmiVector &xBdot,
                     SUNMatrix JB) {
    fJB(t, x.getNVector(), xB.getNVector(), xBdot.getNVector(), JB);
}

void Model_ODE::fJB(realtype t, N_Vector x, N_Vector xB, N_Vector /*xBdot*/,
                    SUNMatrix JB) {
    fJSparse(t, x, J_.get());
    J_.transpose(JB, -1.0, nxtrue_solver);
}

void Model_ODE::fJSparseB(const realtype t, realtype /*cj*/, const AmiVector &x,
                          const AmiVector & /*dx*/, const AmiVector &xB,
                          const AmiVector & /*dxB*/, const AmiVector &xBdot,
                          SUNMatrix JB) {
    fJSparseB(t, x.getNVector(), xB.getNVector(), xBdot.getNVector(), JB);
}

void Model_ODE::fJSparseB(realtype t, N_Vector x, N_Vector xB,
                          N_Vector /*xBdot*/, SUNMatrix JB) {
    fJSparse(t, x, J_.get());
    J_.transpose(JB, -1.0, nxtrue_solver);
}

void Model_ODE::fJDiag(realtype t, N_Vector JDiag, N_Vector x) {
    fJSparse(t, x, J_.get());
    J_.to_diag(JDiag);
}

void Model_ODE::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                     N_Vector xB) {
    N_VConst(0.0, JvB);
    fJSparseB(t, x, xB, nullptr, JB_.get());
    JB_.multiply(JvB, vB);
}

void Model_ODE::fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot) {
    N_VConst(0.0, xBdot);
    fJSparseB(t, x, xB, nullptr, JB_.get());
    JB_.multiply(xBdot, xB);
}

void Model_ODE::fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot) {
    /* initialize with zeros */
    N_VConst(0.0, qBdot);
    fdxdotdp(t, x);

    if (pythonGenerated) {
        /* call multiplication */
        dxdotdp_full.multiply(qBdot, xB, state_.plist, true);
        N_VScale(-1.0, qBdot, qBdot);
    } else {
        /* was matlab generated */
        for (int ip = 0; ip < nplist(); ip++) {
            for (int ix = 0; ix < nxtrue_solver; ix++)
                NV_Ith_S(qBdot, ip * nJ) -= NV_Ith_S(xB, ix) * dxdotdp.at(ix, ip);
            // second order part
            for (int iJ = 1; iJ < nJ; iJ++)
                for (int ix = 0; ix < nxtrue_solver; ix++)
                    NV_Ith_S(qBdot, ip * nJ + iJ) -=
                    NV_Ith_S(xB, ix) * dxdotdp.at(ix + iJ * nxtrue_solver, ip) +
                    NV_Ith_S(xB, ix + iJ * nxtrue_solver) * dxdotdp.at(ix, ip);
        }
    }
}

void Model_ODE::fxBdot_ss(const realtype t, const AmiVector &xB,
                          const AmiVector & /*dx*/, AmiVector &xBdot) {
    fxBdot_ss(t, xB.getNVector(), xBdot.getNVector());
}

void Model_ODE::fxBdot_ss(realtype /*t*/, N_Vector xB, N_Vector xBdot) const {
    /* Right hande side of the adjoint state for steady state computations.
       J is fixed (as x remeins in steady state), so the RHS becomes simple. */
    N_VConst(0.0, xBdot);
    JB_.multiply(xBdot, xB);
}

void Model_ODE::fqBdot_ss(realtype /*t*/, N_Vector xB, N_Vector qBdot) const {
    /* Quadratures when computing adjoints for steady state. The integrand is
       just the adjoint state itself. */
    N_VScale(1.0, xB, qBdot);
}

void Model_ODE::fJSparseB_ss(SUNMatrix JB) {
    /* Just copy the model Jacobian */
    SUNMatCopy(JB_.get(), JB);
}

void Model_ODE::writeSteadystateJB(const realtype t, realtype /*cj*/,
                                   const AmiVector &x, const AmiVector & /*dx*/,
                                   const AmiVector &xB, const AmiVector & /*dxB*/,
                                   const AmiVector &xBdot) {
    /* Get backward Jacobian */
    fJSparseB(t, x.getNVector(), xB.getNVector(), xBdot.getNVector(), JB_.get());
    /* Switch sign, as we integrate forward in time, not backward */
    JB_.scale(-1);
}

void Model_ODE::fsxdot(const realtype t, const AmiVector &x,
                       const AmiVector & /*dx*/, const int ip,
                       const AmiVector &sx, const AmiVector & /*sdx*/,
                       AmiVector &sxdot) {
    fsxdot(t, x.getNVector(), ip, sx.getNVector(), sxdot.getNVector());
}

void Model_ODE::fsxdot(realtype t, N_Vector x, int ip, N_Vector sx,
                       N_Vector sxdot) {

    /* sxdot is just the total derivative d(xdot)dp,
     so we just call dxdotdp and copy the stuff over */
    if (ip == 0) {
        // we only need to call this for the first parameter index will be
        // the same for all remaining
        fdxdotdp(t, x);
        fJSparse(t, x, J_.get());
    }
    if (pythonGenerated) {
        /* copy dxdotdp and the implicit version over */
        // initialize
        N_VConst(0.0, sxdot);
        realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);

        dxdotdp_full.scatter(plist(ip), 1.0, nullptr,
                             gsl::make_span(sxdot_tmp, nx_solver),
                             0, nullptr, 0);

    } else {
        /* copy dxdotdp over */
        N_VScale(1.0, dxdotdp.getNVector(ip), sxdot);
    }
    J_.multiply(sxdot, sx);
}

} // namespace amici
