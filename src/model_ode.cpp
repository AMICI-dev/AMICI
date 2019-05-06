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
    SUNMatZero(J);
    fJ(SM_DATA_D(J), t, N_VGetArrayPointer(x_pos), unscaledParameters.data(),
       fixedParameters.data(), h.data(), w.data(), dwdx.data());
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
    if (wasPythonGenerated()) {
        fJSparse(SM_DATA_S(J), t, N_VGetArrayPointer(x_pos),
                 unscaledParameters.data(), fixedParameters.data(), h.data(),
                 w.data(), dwdx.data());
        fJSparse_colptrs(SM_INDEXPTRS_S(J));
        fJSparse_rowvals(SM_INDEXVALS_S(J));
    } else {
        fJSparse(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(J)), t,
                 N_VGetArrayPointer(x_pos), unscaledParameters.data(),
                 fixedParameters.data(), h.data(), w.data(), dwdx.data());
    }
}

void Model_ODE::fJv(const realtype t, const AmiVector &x,
                    const AmiVector & /*dx*/, const AmiVector & /*xdot*/,
                    const AmiVector &v, AmiVector &Jv, const realtype /*cj*/) {
    fJv(v.getNVector(), Jv.getNVector(), t, x.getNVector());
}

void Model_ODE::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x) {
    N_VConst(0.0, Jv);
    fJSparse(t, x, J.get());
    J.multiply(Jv, v);
}

void Model_ODE::froot(const realtype t, const AmiVector &x,
                      const AmiVector & /*dx*/, gsl::span<realtype> root) {
    froot(t, x.getNVector(), root);
}

void Model_ODE::froot(realtype t, N_Vector x, gsl::span<realtype> root) {
    auto x_pos = computeX_pos(x);
    std::fill(root.begin(), root.end(), 0.0);
    froot(root.data(), t, N_VGetArrayPointer(x_pos), unscaledParameters.data(),
          fixedParameters.data(), h.data());
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
          unscaledParameters.data(), fixedParameters.data(), h.data(),
          w.data());
}

void Model_ODE::fJDiag(const realtype t, AmiVector &JDiag,
                       const realtype /*cj*/, const AmiVector &x,
                       const AmiVector & /*dx*/) {
    fJDiag(t, JDiag.getNVector(), x.getNVector());
    if (checkFinite(JDiag.getVector(), "Jacobian") != AMICI_SUCCESS)
        throw AmiException("Evaluation of fJDiag failed!");
}

void Model_ODE::fdxdotdw(const realtype t, const N_Vector x) {
    dxdotdw.reset();
    auto x_pos = computeX_pos(x);
    fdxdotdw(dxdotdw.data(), t, N_VGetArrayPointer(x_pos),
             unscaledParameters.data(), fixedParameters.data(), h.data(),
             w.data());
    fdxdotdw_colptrs(dxdotdw.indexptrs());
    fdxdotdw_rowvals(dxdotdw.indexvals());
}

void Model_ODE::fdxdotdp(const realtype t, const N_Vector x) {
    auto x_pos = computeX_pos(x);
    fdwdp(t, N_VGetArrayPointer(x));
    if (wasPythonGenerated()) {
        // python generated
        fdxdotdw(t, x);
        for (int ip = 0; ip < nplist(); ip++) {
            N_VConst(0.0, dxdotdp.getNVector(ip));
            fdxdotdp(dxdotdp.data(ip), t, N_VGetArrayPointer(x_pos),
                     unscaledParameters.data(), fixedParameters.data(),
                     h.data(), plist_[ip], w.data());
            if (nw > 0)
                dxdotdw.multiply(
                    gsl::span<realtype>(dxdotdp.data(ip), nx_solver),
                    gsl::span<const realtype>(&dwdp.at(nw * ip), nw));
        }
    } else {
        // matlab generated
        for (int ip = 0; ip < nplist(); ip++) {
            N_VConst(0.0, dxdotdp.getNVector(ip));
            fdxdotdp(dxdotdp.data(ip), t, N_VGetArrayPointer(x_pos),
                     unscaledParameters.data(), fixedParameters.data(),
                     h.data(), plist_[ip], w.data(), dwdp.data());
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

void Model_ODE::fJB(realtype * /*JB*/, const realtype /*t*/,
                    const realtype * /*x*/, const realtype * /*p*/,
                    const realtype * /*k*/, const realtype * /*h*/,
                    const realtype * /*xB*/, const realtype * /*w*/,
                    const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s is not "
                       "implemented for this model!",
                       __func__);
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

void Model_ODE::fJSparseB(SUNMatrixContent_Sparse /*JSparseB*/,
                          const realtype /*t*/, const realtype * /*x*/,
                          const realtype * /*p*/, const realtype * /*k*/,
                          const realtype * /*h*/, const realtype * /*xB*/,
                          const realtype * /*w*/, const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparseB(realtype * /*JSparseB*/, const realtype /*t*/,
                          const realtype * /*x*/, const realtype * /*p*/,
                          const realtype * /*k*/, const realtype * /*h*/,
                          const realtype * /*xB*/, const realtype * /*w*/,
                          const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparseB_colptrs(sunindextype * /*indexptrs*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJSparseB_rowvals(sunindextype * /*indexvals*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__); // not implemented
}

void Model_ODE::fJDiag(realtype * /*JDiag*/, const realtype /*t*/,
                       const realtype * /*x*/, const realtype * /*p*/,
                       const realtype * /*k*/, const realtype * /*h*/,
                       const realtype * /*w*/, const realtype * /*dwdx*/) {
    throw AmiException("Requested functionality is not supported as %s is not "
                       "implemented for this model!",
                       __func__);
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

void Model_ODE::fdxdotdp(realtype * /*dxdotdp*/, const realtype /*t*/,
                         const realtype * /*x*/, const realtype * /*p*/,
                         const realtype * /*k*/, const realtype * /*h*/,
                         const int /*ip*/, const realtype * /*w*/) {
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

void Model_ODE::fJB(realtype t, N_Vector x, N_Vector xB, N_Vector /*xBdot*/,
                    SUNMatrix JB) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointer(x_pos));
    SUNMatZero(JB);
    fJB(SM_DATA_D(JB), t, N_VGetArrayPointer(x_pos), unscaledParameters.data(),
        fixedParameters.data(), h.data(), N_VGetArrayPointer(xB), w.data(),
        dwdx.data());
}

void Model_ODE::fJSparseB(realtype t, N_Vector x, N_Vector xB,
                          N_Vector /*xBdot*/, SUNMatrix JB) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointer(x_pos));
    SUNMatZero(JB);
    if (wasPythonGenerated()) {
        fJSparseB(SM_DATA_S(JB), t, N_VGetArrayPointer(x_pos),
                  unscaledParameters.data(), fixedParameters.data(), h.data(),
                  N_VGetArrayPointer(xB), w.data(), dwdx.data());
        fJSparseB_colptrs(SM_INDEXPTRS_S(JB));
        fJSparseB_rowvals(SM_INDEXVALS_S(JB));
    } else {
        fJSparseB(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(JB)), t,
                  N_VGetArrayPointer(x_pos), unscaledParameters.data(),
                  fixedParameters.data(), h.data(), N_VGetArrayPointer(xB),
                  w.data(), dwdx.data());
    }
}

void Model_ODE::fJDiag(realtype t, N_Vector JDiag, N_Vector x) {
    auto x_pos = computeX_pos(x);
    fdwdx(t, N_VGetArrayPointer(x_pos));
    N_VConst(0.0, JDiag);
    fJDiag(N_VGetArrayPointer(JDiag), t, N_VGetArrayPointer(x_pos),
           unscaledParameters.data(), fixedParameters.data(), h.data(),
           w.data(), dwdx.data());
}

void Model_ODE::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                     N_Vector xB) {
    N_VConst(0.0, JvB);
    fJSparseB(t, x, xB, nullptr, J.get());
    J.multiply(JvB, vB);
}

void Model_ODE::fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot) {
    N_VConst(0.0, xBdot);
    fJSparseB(t, x, xB, nullptr, J.get());
    J.multiply(xBdot, xB);
}

void Model_ODE::fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot) {
    N_VConst(0.0, qBdot);
    fdxdotdp(t, x);
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

void Model_ODE::fsxdot(const realtype t, const AmiVector &x,
                       const AmiVector & /*dx*/, const int ip,
                       const AmiVector &sx, const AmiVector & /*sdx*/,
                       AmiVector &sxdot) {
    fsxdot(t, x.getNVector(), ip, sx.getNVector(), sxdot.getNVector());
}

void Model_ODE::fsxdot(realtype t, N_Vector x, int ip, N_Vector sx,
                       N_Vector sxdot) {
    if (ip == 0) {
        // we only need to call this for the first parameter index will be
        // the same for all remaining
        fdxdotdp(t, x);
        fJSparse(t, x, J.get());
    }
    N_VScale(1.0, dxdotdp.getNVector(ip), sxdot);
    J.multiply(sxdot, sx);
}

} // namespace amici
