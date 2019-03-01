#include "amici/solver_cvodes.h"
#include "amici/model_ode.h"

namespace amici {

    void Model_ODE::fJ(realtype t, realtype  /*cj*/, AmiVector *x, AmiVector * /*dx*/,
                          AmiVector *xdot, SUNMatrix J) {
        fJ(t, x->getNVector(), xdot->getNVector(), J);

    }

    /** implementation of fJ at the N_Vector level, this function provides an
     *interface to the model specific routines for the solver implementation
     *aswell as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     **/
    void Model_ODE::fJ(realtype t, N_Vector x, N_Vector  /*xdot*/, SUNMatrix J) {
        auto x_pos = computeX_pos(x);
        fdwdx(t, N_VGetArrayPointer(x_pos));
        SUNMatZero(J);
        fJ(SM_DATA_D(J), t, N_VGetArrayPointer(x_pos),
           unscaledParameters.data(), fixedParameters.data(), h.data(),
           w.data(), dwdx.data());
    }

    void Model_ODE::fJSparse(realtype t, realtype  /*cj*/, AmiVector *x,
                             AmiVector * /*dx*/, AmiVector * /*xdot*/, SUNMatrix J) {
        fJSparse(t, x->getNVector(), J);
    }

    /** implementation of fJSparse at the N_Vector level, this function provides
     * an interface to the model specific routines for the solver implementation
     * aswell as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param J Matrix to which the Jacobian will be written
     */
    void Model_ODE::fJSparse(realtype t, N_Vector x, SUNMatrix J) {
        auto x_pos = computeX_pos(x);
        fdwdx(t, N_VGetArrayPointer(x_pos));
        SUNMatZero(J);
        if (wasPythonGenerated()) {
            fJSparse(SM_DATA_S(J), t, N_VGetArrayPointer(x_pos),
                     unscaledParameters.data(), fixedParameters.data(),
                     h.data(), w.data(), dwdx.data());
            fJSparse_colptrs(SM_INDEXPTRS_S(J));
            fJSparse_rowvals(SM_INDEXVALS_S(J));
        } else {
            fJSparse(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(J)), t,
                     N_VGetArrayPointer(x_pos), unscaledParameters.data(),
                     fixedParameters.data(), h.data(), w.data(), dwdx.data());
        }
    }

    void Model_ODE::fJv(realtype t, AmiVector *x, AmiVector * /*dx*/,
                        AmiVector * /*xdot*/, AmiVector *v, AmiVector *Jv,
                        realtype /*cj*/) {
        fJv(v->getNVector(), Jv->getNVector(), t, x->getNVector());
    }

    /** implementation of fJv at the N_Vector level.
     * @param t timepoint
     * @param x Vector with the states
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     * written
     **/
    void Model_ODE::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x) {
        N_VConst(0.0, Jv);
        fJSparse(t, x, J.get());
        J.multiply(Jv, v);
    }

    void Model_ODE::froot(realtype t, AmiVector *x, AmiVector * /*dx*/, realtype *root){
        froot(t,x->getNVector(),root);
    }

    /** implementation of froot at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation aswell as the AmiVector
     * level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param root array with root function values
     */
    void Model_ODE::froot(realtype t, N_Vector x, realtype *root) {
        auto x_pos = computeX_pos(x);
        memset(root,0,sizeof(realtype)*ne);
        froot(root,t,N_VGetArrayPointer(x_pos),unscaledParameters.data(),fixedParameters.data(),h.data());
    }

    void Model_ODE::fxdot(realtype t, AmiVector *x, AmiVector * /*dx*/, AmiVector *xdot) {
        fxdot(t,x->getNVector(),xdot->getNVector());
    }

    /** implementation of fxdot at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation aswell as the AmiVector
     * level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     */
    void Model_ODE::fxdot(realtype t, N_Vector x, N_Vector xdot) {
        auto x_pos = computeX_pos(x);
        fw(t,N_VGetArrayPointer(x_pos));
        N_VConst(0.0,xdot);
        fxdot(N_VGetArrayPointer(xdot),t,N_VGetArrayPointer(x_pos),unscaledParameters.data(),fixedParameters.data(),h.data(),
                     w.data());
    }

    /** diagonalized Jacobian (for preconditioning)
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @return status flag indicating successful execution
     **/
    void Model_ODE::fJDiag(realtype t, AmiVector *JDiag, realtype  /*cj*/, AmiVector *x,
                          AmiVector * /*dx*/) {
        fJDiag(t, JDiag->getNVector(), x->getNVector());
        if(checkFinite(nx_solver,JDiag->data(),"Jacobian") != AMICI_SUCCESS)
            throw AmiException("Evaluation of fJDiag failed!");
    }

    /** Sensitivity of dx/dt wrt model parameters w
     * @param t timepoint
     * @param x Vector with the states
     * @return status flag indicating successful execution
     */
    void Model_ODE::fdxdotdw(const realtype t, const N_Vector x) {
        dxdotdw.reset();
        auto x_pos = computeX_pos(x);
        fdxdotdw(dxdotdw.data(), t, N_VGetArrayPointer(x_pos),
                 unscaledParameters.data(), fixedParameters.data(),
                 h.data(), w.data());
        fdxdotdw_colptrs(dxdotdw.indexptrs());
        fdxdotdw_rowvals(dxdotdw.indexvals());
    }

    /** Sensitivity of dx/dt wrt model parameters p
     * @param t timepoint
     * @param x Vector with the states
     * @return status flag indicating successful execution
     */
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
                    dxdotdw.multiply(dxdotdp.data(ip), &dwdp.at(nw * ip));
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

    std::unique_ptr<Solver> Model_ODE::getSolver() {
        return std::unique_ptr<Solver>(new amici::CVodeSolver());
    }

    /** implementation of fJB at the N_Vector level, this function provides an
     *interface to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     **/
    void Model_ODE::fJB(realtype t, N_Vector x, N_Vector xB, N_Vector  /*xBdot*/,
                        SUNMatrix JB) {
        auto x_pos = computeX_pos(x);
        fdwdx(t, N_VGetArrayPointer(x_pos));
        SUNMatZero(JB);
        fJB(SM_DATA_D(JB), t, N_VGetArrayPointer(x_pos),
            unscaledParameters.data(), fixedParameters.data(), h.data(),
            N_VGetArrayPointer(xB), w.data(), dwdx.data());
    }

    /** implementation of fJSparseB at the N_Vector level, this function
     * provides an interface to the model specific routines for the solver
     * implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     */
    void Model_ODE::fJSparseB(realtype t, N_Vector x, N_Vector xB,
                              N_Vector /*xBdot*/, SUNMatrix JB) {
        auto x_pos = computeX_pos(x);
        fdwdx(t, N_VGetArrayPointer(x_pos));
        SUNMatZero(JB);
        if (wasPythonGenerated()) {
            fJSparseB(SM_DATA_S(JB), t, N_VGetArrayPointer(x_pos),
                      unscaledParameters.data(), fixedParameters.data(),
                      h.data(), N_VGetArrayPointer(xB), w.data(), dwdx.data());
            fJSparseB_colptrs(SM_INDEXPTRS_S(JB));
            fJSparseB_rowvals(SM_INDEXVALS_S(JB));
        } else {
            fJSparseB(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(JB)), t,
                      N_VGetArrayPointer(x_pos), unscaledParameters.data(),
                      fixedParameters.data(), h.data(), N_VGetArrayPointer(xB),
                      w.data(), dwdx.data());
        }
    }

    /** implementation of fJDiag at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param x Vector with the states
     **/
    void Model_ODE::fJDiag(realtype t, N_Vector JDiag, N_Vector x) {
        auto x_pos = computeX_pos(x);
        fdwdx(t,N_VGetArrayPointer(x_pos));
        N_VConst(0.0,JDiag);
        fJDiag(N_VGetArrayPointer(JDiag),t,N_VGetArrayPointer(x_pos),unscaledParameters.data(),fixedParameters.data(),h.data(),
                      w.data(),dwdx.data());
    }

    /** implementation of fJvB at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be
     *written
     **/
    void Model_ODE::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                         N_Vector xB) {
        N_VConst(0.0, JvB);
        fJSparseB(t, x, xB, nullptr, J.get());
        J.multiply(JvB, vB);
    }

    /** implementation of fxBdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     */
    void Model_ODE::fxBdot(realtype t, N_Vector x, N_Vector xB,
                           N_Vector xBdot) {
        N_VConst(0.0, xBdot);
        fJSparseB(t, x, xB, nullptr, J.get());
        J.multiply(xBdot, xB);
    }

    /** implementation of fqBdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void Model_ODE::fqBdot(realtype t, N_Vector x, N_Vector xB,
                           N_Vector qBdot) {
        N_VConst(0.0, qBdot);
        fdxdotdp(t, x);
        for (int ip = 0; ip < nplist(); ip++) {
            for (int ix = 0; ix < nxtrue_solver; ix++)
                NV_Ith_S(qBdot, ip * nJ) -=
                NV_Ith_S(xB, ix) * dxdotdp.at(ix, ip);
            // second order part
            for (int iJ = 1; iJ < nJ; iJ++)
                for (int ix = 0; ix < nxtrue_solver; ix++)
                    NV_Ith_S(qBdot, ip * nJ + iJ) -=
                        NV_Ith_S(xB, ix) *
                            dxdotdp.at(ix + iJ * nxtrue_solver, ip) +
                        NV_Ith_S(xB, ix + iJ * nxtrue_solver) *
                            dxdotdp.at(ix, ip);
        }
    }

    void Model_ODE::fsxdot(realtype t, AmiVector *x, AmiVector * /*dx*/, int ip,
                           AmiVector *sx, AmiVector * /*sdx*/, AmiVector *sxdot) {
        fsxdot(t,x->getNVector(), ip, sx->getNVector(), sxdot->getNVector());
    }

    /** implementation of fsxdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     */
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
}
