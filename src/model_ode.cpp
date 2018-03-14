#include "amici/solver_cvodes.h"
#include "amici/model_ode.h"

namespace amici {
    
    void Model_ODE::fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                          AmiVector *xdot, DlsMat J) {
        fJ(t, x->getNVector(), xdot->getNVector(), J);
        
    }
    
    /** implementation of fJ at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation aswell as the AmiVector
     * level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     **/
    void Model_ODE::fJ( realtype t, N_Vector x, N_Vector xdot, DlsMat J) {
        fdwdx(t,x);
        SetToZero(J);
        fJ(J->data,t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
           w.data(),dwdx.data());
    }
    
    void Model_ODE::fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                        AmiVector *xdot, SlsMat J){
        fJSparse(t,x->getNVector(),J);
    }
    
    /** implementation of fJSparse at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation aswell as the AmiVector
     * level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param J Matrix to which the Jacobian will be written
     */
    void Model_ODE::fJSparse(realtype t, N_Vector x, SlsMat J) {
        fdwdx(t,x);
        SparseSetMatToZero(J);
        fJSparse(J,t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                 w.data(),dwdx.data());
    }
    
    void Model_ODE::fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                        AmiVector *v, AmiVector *Jv, realtype cj){
        fJv(v->getNVector(),Jv->getNVector(),t,x->getNVector());
    }
    
    /** implementation of fJv at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation aswell as the AmiVector
     * level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     * written
     **/
    void Model_ODE::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x) {
        fdwdx(t,x);
        N_VConst(0.0,Jv);
        fJv(N_VGetArrayPointer(Jv),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
            N_VGetArrayPointer(v),w.data(),dwdx.data());
    }
    
    void Model_ODE::froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root){
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
        memset(root,0.0,sizeof(realtype)*ne);
        froot(root,t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data());
    }
    
    void Model_ODE::fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot) {
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
        fw(t,x);
        N_VConst(0.0,xdot);
        fxdot(N_VGetArrayPointer(xdot),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
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
    void Model_ODE::fJDiag(realtype t, AmiVector *JDiag, realtype cj, AmiVector *x,
                          AmiVector *dx) {
        fJDiag(t, JDiag->getNVector(), x->getNVector());
        if(!checkFinite(nx,JDiag->data(),"Jacobian"))
            throw AmiException("Evaluation of fJDiag failed!");
    }
    
    /** Sensitivity of dx/dt wrt model parameters p
     * @param t timepoint
     * @param x Vector with the states
     * @return status flag indicating successful execution
     */
    void Model_ODE::fdxdotdp(const realtype t, const N_Vector x) {
        std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
        fdwdp(t,x);
        for(int ip = 0; ip < nplist(); ip++)
            fdxdotdp(&dxdotdp.at(nx*ip),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                     plist_[ip],w.data(),dwdp.data());
    }
    
    std::unique_ptr<Solver> Model_ODE::getSolver() {
        return std::unique_ptr<Solver>(new amici::CVodeSolver());
    }
    
    /** implementation of fJB at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     **/
    void Model_ODE::fJB( realtype t, N_Vector x, N_Vector xB,
                         N_Vector xBdot, DlsMat JB) {
        fdwdx(t,x);
        SetToZero(JB);
        fJB(JB->data,t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                   N_VGetArrayPointer(xB),w.data(),dwdx.data());
    }
    
    /** implementation of fJSparseB at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     */
    void Model_ODE::fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB) {
        fdwdx(t,x);
        SparseSetMatToZero(JB);
        fJSparseB(JB,t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                         N_VGetArrayPointer(xB),w.data(),dwdx.data());
    }
    
    /** implementation of fJDiag at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param x Vector with the states
     **/
    void Model_ODE::fJDiag(realtype t, N_Vector JDiag, N_Vector x) {
        fdwdx(t,x);
        N_VConst(0.0,JDiag);
        fJDiag(N_VGetArrayPointer(JDiag),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                      w.data(),dwdx.data());
    }
    
    /** implementation of fJvB at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be
     *written
     **/
    void Model_ODE::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB) {
        fdwdx(t,x);
        N_VConst(0.0,JvB);
        fJvB(N_VGetArrayPointer(JvB),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                    N_VGetArrayPointer(xB),N_VGetArrayPointer(vB),w.data(),dwdx.data());
    }
    
    /** implementation of fxBdot at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     */
    void Model_ODE::fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot) {
        fdwdx(t,x);
        N_VConst(0.0,xBdot);
        fxBdot(N_VGetArrayPointer(xBdot),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                      N_VGetArrayPointer(xB),w.data(),dwdx.data());
    }
    
    /** implementation of fqBdot at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void Model_ODE::fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot) {
        fdwdp(t,x);
        N_VConst(0.0,qBdot);
        realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
        for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
            fqBdot(&qBdot_tmp[ip*nJ],plist_[ip],t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                          N_VGetArrayPointer(xB),w.data(),dwdp.data());
    }
    
    /** implementation of fsxdot at the N_Vector level, this function provides an interface
     * to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     */
    void Model_ODE::fsxdot(realtype t, N_Vector x, int ip,
                            N_Vector sx, N_Vector sxdot) {
        if(ip == 0) { // we only need to call this for the first parameter index will be the same for all remaining
            fdxdotdp(t, x);
            fJSparse(t, x, J);
        }
        N_VConst(0.0,sxdot);
        fsxdot(N_VGetArrayPointer(sxdot),t,N_VGetArrayPointer(x), unscaledParameters.data(),fixedParameters.data(),h.data(),
                      plist_[ip],N_VGetArrayPointer(sx),
                      w.data(),dwdx.data(),J->data,&dxdotdp.at(ip*nx));
    }
}
