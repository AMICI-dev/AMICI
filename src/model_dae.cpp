#include "amici/model_dae.h"
#include "amici/solver_idas.h"

namespace amici {

void Model_DAE::fJ(const realtype t, const realtype cj, const AmiVector &x,
                   const AmiVector &dx, const AmiVector &xdot, SUNMatrix J) {
    fJ(t, cj, x.getNVector(), dx.getNVector(), xdot.getNVector(), J);
}

void Model_DAE::fJ(realtype t, realtype cj, const_N_Vector x, const_N_Vector dx,
                   const_N_Vector /*xdot*/, SUNMatrix J) {
    fJSparse(t, cj, x, dx, derived_state_.J_.get());
    derived_state_.J_.refresh();
    auto JDense = SUNMatrixWrapper(J);
    derived_state_.J_.to_dense(JDense);
}

void Model_DAE::fJSparse(const realtype t, const realtype cj,
                         const AmiVector &x, const AmiVector &dx,
                         const AmiVector & /*xdot*/, SUNMatrix J) {
    fJSparse(t, cj, x.getNVector(), dx.getNVector(), J);
}

void Model_DAE::fJSparse(realtype t, realtype cj, const_N_Vector x,
                         const_N_Vector dx, SUNMatrix J) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointerConst(x_pos));
    SUNMatZero(J);
    fJSparse(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(J)), t,
             N_VGetArrayPointerConst(x_pos),
             state_.unscaledParameters.data(),
             state_.fixedParameters.data(), state_.h.data(), cj,
             N_VGetArrayPointerConst(dx),
             derived_state_.w_.data(), derived_state_.dwdx_.data());
}

void Model_DAE::fJv(const realtype t, const AmiVector &x, const AmiVector &dx,
                    const AmiVector & /*xdot*/, const AmiVector &v,
                    AmiVector &Jv, const realtype cj) {
    fJv(t, x.getNVector(), dx.getNVector(), v.getNVector(), Jv.getNVector(),
        cj);
}

void Model_DAE::fJv(realtype t, const_N_Vector x, const_N_Vector dx,
                    const_N_Vector v, N_Vector Jv, realtype cj) {
    N_VConst(0.0, Jv);
    fJSparse(t, cj, x, dx, derived_state_.J_.get());
    derived_state_.J_.refresh();
    derived_state_.J_.multiply(Jv, v);
}

void Model_DAE::froot(const realtype t, const AmiVector &x, const AmiVector &dx,
                      gsl::span<realtype> root) {
    froot(t, x.getNVector(), dx.getNVector(), root);
}

void Model_DAE::froot(realtype t, const_N_Vector x, const_N_Vector dx,
                      gsl::span<realtype> root) {
    std::fill(root.begin(), root.end(), 0.0);
    auto x_pos = computeX_pos(x);
    froot(root.data(), t, N_VGetArrayPointerConst(x_pos),
          state_.unscaledParameters.data(), state_.fixedParameters.data(),
          state_.h.data(), N_VGetArrayPointerConst(dx));
}

void Model_DAE::fxdot(const realtype t, const AmiVector &x, const AmiVector &dx,
                      AmiVector &xdot) {
    fxdot(t, x.getNVector(), dx.getNVector(), xdot.getNVector());
}

void Model_DAE::fxdot(realtype t, const_N_Vector x, const_N_Vector dx,
                      N_Vector xdot) {
    auto x_pos = computeX_pos(x);
    fw(t, N_VGetArrayPointerConst(x));
    N_VConst(0.0, xdot);
    fxdot(N_VGetArrayPointer(xdot), t,
          N_VGetArrayPointerConst(x_pos),
          state_.unscaledParameters.data(), state_.fixedParameters.data(),
          state_.h.data(), N_VGetArrayPointerConst(dx),
          derived_state_.w_.data());
}

void Model_DAE::fJDiag(const realtype t, AmiVector &JDiag,
                       const realtype /*cj*/, const AmiVector &x,
                       const AmiVector &dx) {
    fJSparse(t, 0.0, x.getNVector(), dx.getNVector(), derived_state_.J_.get());
    derived_state_.J_.refresh();
    derived_state_.J_.to_diag(JDiag.getNVector());
    if (!checkFinite(JDiag.getVector(), "Jacobian"))
        throw AmiException("Evaluation of fJDiag failed!");
}

void Model_DAE::fdxdotdp(const realtype t, const const_N_Vector x,
                         const const_N_Vector dx) {
    auto x_pos = computeX_pos(x);

    if (pythonGenerated) {
        // python generated, not yet implemented for DAEs
        throw AmiException("Wrapping of DAEs is not yet implemented from Python");
    } else {
        // matlab generated
        fdwdp(t, N_VGetArrayPointerConst(x_pos));

        for (int ip = 0; ip < nplist(); ip++) {
            N_VConst(0.0, derived_state_.dxdotdp.getNVector(ip));
            fdxdotdp(derived_state_.dxdotdp.data(ip), t,
                     N_VGetArrayPointerConst(x_pos),
                     state_.unscaledParameters.data(),
                     state_.fixedParameters.data(), state_.h.data(), plist(ip),
                     N_VGetArrayPointerConst(dx), derived_state_.w_.data(),
                     derived_state_.dwdp_.data());
        }
    }
}

void Model_DAE::fM(realtype t, const_N_Vector x) {
    derived_state_.M_.zero();
    auto x_pos = computeX_pos(x);
    fM(derived_state_.M_.data(), t, N_VGetArrayPointerConst(x_pos),
       state_.unscaledParameters.data(),
       state_.fixedParameters.data());
}

std::unique_ptr<Solver> Model_DAE::getSolver() {
    return std::unique_ptr<Solver>(new amici::IDASolver());
}

void Model_DAE::froot(realtype * /*root*/, const realtype /*t*/,
                      const realtype * /*x*/, const double * /*p*/, const double * /*k*/,
                      const realtype * /*h*/, const realtype * /*dx*/) {
    throw AmiException("Requested functionality is not supported as %s is not "
                       "implemented for this model!",
                       __func__); // not implemented
}

void Model_DAE::fdxdotdp(realtype * /*dxdotdp*/, const realtype /*t*/,
                         const realtype * /*x*/, const realtype * /*p*/,
                         const realtype * /*k*/, const realtype * /*h*/,
                         const int /*ip*/, const realtype * /*dx*/,
                         const realtype * /*w*/, const realtype * /*dwdp*/) {
    throw AmiException("Requested functionality is not supported as %s is not "
                       "implemented for this model!",
                       __func__);
}

void Model_DAE::fM(realtype */*M*/, const realtype /*t*/, const realtype */*x*/,
                   const realtype */*p*/, const realtype */*k*/){}

void Model_DAE::fJB(const realtype t, realtype cj, const AmiVector &x,
                    const AmiVector &dx, const AmiVector &xB,
                     const AmiVector &/*dxB*/, const AmiVector & /*xBdot*/,
                     SUNMatrix JB) {
    fJB(t, cj, x.getNVector(), dx.getNVector(), xB.getNVector(),
        dx.getNVector(), JB);
}


void Model_DAE::fJB(realtype t, realtype cj, const_N_Vector x,
                    const_N_Vector dx, const_N_Vector /*xB*/,
                    const_N_Vector /*dxB*/, SUNMatrix JB) {
    fJSparse(t, cj, x, dx, derived_state_.J_.get());
    derived_state_.J_.refresh();
    auto JBDense = SUNMatrixWrapper(JB);
    derived_state_.J_.transpose(JBDense, -1.0, nxtrue_solver);
}

void Model_DAE::fJSparseB(const realtype t, realtype cj, const AmiVector &x,
                          const AmiVector &dx, const AmiVector &xB,
                          const AmiVector &dxB, const AmiVector & /*xBdot*/,
                          SUNMatrix JB) {
    fJSparseB(t, cj, x.getNVector(), dx.getNVector(), xB.getNVector(), dxB.getNVector(), JB);
}

void Model_DAE::fJSparseB(realtype t, realtype cj, const_N_Vector x,
                          const_N_Vector dx,
                          const_N_Vector /*xB*/, const_N_Vector /*dxB*/,
                          SUNMatrix JB) {
    fJSparse(t, cj, x, dx, derived_state_.J_.get());
    derived_state_.J_.refresh();
    auto JSparseB = SUNMatrixWrapper(JB);
    derived_state_.J_.transpose(JSparseB, -1.0, nxtrue_solver);
}

void Model_DAE::fJvB(realtype t, const_N_Vector x, const_N_Vector dx,
                     const_N_Vector xB, const_N_Vector dxB, const_N_Vector vB,
                     N_Vector JvB, realtype cj) {
    N_VConst(0.0, JvB);
    fJSparseB(t, cj, x, dx, xB, dxB, derived_state_.JB_.get());
    derived_state_.JB_.refresh();
    derived_state_.JB_.multiply(JvB, vB);
}

void Model_DAE::fxBdot(realtype t, const_N_Vector x, const_N_Vector dx,
                       const_N_Vector xB,
                       const_N_Vector dxB, N_Vector xBdot) {
    N_VConst(0.0, xBdot);
    fJSparseB(t, 1.0, x, dx, xB, dxB, derived_state_.JB_.get());
    derived_state_.JB_.refresh();
    fM(t, x);
    derived_state_.JB_.multiply(xBdot, xB);
}

void Model_DAE::fqBdot(realtype t, const_N_Vector x, const_N_Vector dx,
                       const_N_Vector xB, const_N_Vector /*dxB*/,
                       N_Vector qBdot) {
    N_VConst(0.0, qBdot);
    fdxdotdp(t, x, dx);
    for (int ip = 0; ip < nplist(); ip++) {
        for (int ix = 0; ix < nxtrue_solver; ix++)
            NV_Ith_S(qBdot, ip * nJ) -= NV_Ith_S(xB, ix)
                    * derived_state_.dxdotdp.at(ix, ip);
        // second order part
        for (int iJ = 1; iJ < nJ; iJ++)
            for (int ix = 0; ix < nxtrue_solver; ix++)
                NV_Ith_S(qBdot, ip * nJ + iJ) -=
                    NV_Ith_S(xB, ix)
                        * derived_state_.dxdotdp.at(ix + iJ * nxtrue_solver, ip)
                        + NV_Ith_S(xB, ix + iJ * nxtrue_solver)
                        * derived_state_.dxdotdp.at(ix, ip);
    }
}

void Model_DAE::fxBdot_ss(const realtype t, const AmiVector &xB,
                          const AmiVector &dxB, AmiVector &xBdot) {
    fxBdot_ss(t, xB.getNVector(), dxB.getNVector(), xBdot.getNVector());
}

void Model_DAE::fxBdot_ss(realtype /*t*/, const_N_Vector xB, const_N_Vector /*dxB*/,
                          N_Vector xBdot) const {
    /* Right hand side of the adjoint state for steady state computations.
     J is fixed (as x remains in steady state), so the RHS becomes simple. */
    N_VConst(0.0, xBdot);
    derived_state_.JB_.multiply(xBdot, xB);
    /* Mind the minus sign... */
    N_VScale(-1.0, xBdot, xBdot);
}

void Model_DAE::fqBdot_ss(realtype /*t*/, const_N_Vector xB,
                          const_N_Vector /*dxB*/, N_Vector qBdot) const {
    /* Quadratures when computing adjoints for steady state. The integrand is
     just the adjoint state itself. */
    N_VScale(1.0, const_cast<N_Vector>(xB), qBdot);
}

void Model_DAE::fJSparseB_ss(SUNMatrix JB) {
    /* Just pass the model Jacobian on to JB */
    SUNMatCopy(derived_state_.JB_.get(), JB);
    derived_state_.JB_.refresh();
}

void Model_DAE::writeSteadystateJB(const realtype t, realtype cj,
                                   const AmiVector &x, const AmiVector & dx,
                                   const AmiVector &xB, const AmiVector & dxB,
                                   const AmiVector &/*xBdot*/) {
    /* Get backward Jacobian */
    fJSparseB(t, cj, x.getNVector(), dx.getNVector(), xB.getNVector(),
              dxB.getNVector(), derived_state_.JB_.get());
    derived_state_.JB_.refresh();
    /* Switch sign, as we integrate forward in time, not backward */
    derived_state_.JB_.scale(-1);
}

void Model_DAE::fsxdot(const realtype t, const AmiVector &x,
                       const AmiVector &dx, const int ip, const AmiVector &sx,
                       const AmiVector &sdx, AmiVector &sxdot) {
    fsxdot(t, x.getNVector(), dx.getNVector(), ip, sx.getNVector(),
           sdx.getNVector(), sxdot.getNVector());
}

void Model_DAE::fsxdot(realtype t, const_N_Vector x, const_N_Vector dx, int ip,
                       const_N_Vector sx, const_N_Vector sdx, N_Vector sxdot) {
    if (ip == 0) {
        // we only need to call this for the first parameter index will be
        // the same for all remaining
        fM(t, x);
        fdxdotdp(t, x, dx);
        fJSparse(t, 0.0, x, dx, derived_state_.J_.get());
        derived_state_.J_.refresh();
    }

    if (pythonGenerated) {
        // python generated, not yet implemented for DAEs
        throw AmiException("Wrapping of DAEs is not yet implemented from Python");
    } else {
        /* copy dxdotdp over */
        N_VScale(1.0, derived_state_.dxdotdp.getNVector(ip), sxdot);
    }

    derived_state_.J_.multiply(sxdot, sx);
    derived_state_.M_.multiply(sxdot, sdx, -1.0);
}

} // namespace amici
