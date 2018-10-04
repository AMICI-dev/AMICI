#include "amici/solver_idas.h"
#include "amici/model_dae.h"

namespace amici {
    
    void Model_DAE::fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                       AmiVector *xdot, DlsMat J) {
        fJ(t,cj, x->getNVector(), dx->getNVector(), xdot->getNVector(), J);
    }
    
    /** Jacobian of xdot with respect to states x
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written

     **/
    void Model_DAE::fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
                      N_Vector xdot, DlsMat J) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        SetToZero(J);
        fJ(J->data,t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                  cj,N_VGetArrayPointer(dx),w.data(),dwdx.data());
    }

     void Model_DAE::fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                              AmiVector *xdot, SlsMat J){
        fJSparse(t,cj,x->getNVector(),dx->getNVector(),J);
    }
    
    /** J in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param J Matrix to which the Jacobian will be written
     */
    void Model_DAE::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, SlsMat J) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        SparseSetMatToZero(J);
        fJSparse(J,t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                        cj,N_VGetArrayPointer(dx),w.data(),dwdx.data());
    }
    
    void Model_DAE::fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                        AmiVector *v, AmiVector *Jv, realtype cj){
        fJv(t,x->getNVector(),dx->getNVector(),v->getNVector(),
                    Jv->getNVector(),cj);
    }
    
    /** Matrix vector product of J with a vector v (for iterative solvers)
     * @param t timepoint @type realtype
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     *written
     **/
    void Model_DAE::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector v, N_Vector Jv,
                       realtype cj) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        N_VConst(0.0,Jv);
        fJv(N_VGetArrayPointer(Jv),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                   cj,N_VGetArrayPointer(dx),N_VGetArrayPointer(v),w.data(),dwdx.data());
    }
    
     void Model_DAE::froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root){
        froot(t,x->getNVector(),dx->getNVector(),root);
    }
    
    /** Event trigger function for events
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param root array with root function values
     */
    void Model_DAE::froot(realtype t, N_Vector x, N_Vector dx, realtype *root) {
        memset(root,0.0,sizeof(realtype)*ne);
        computeX_pos(x);
        froot(root,t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                     N_VGetArrayPointer(dx));
    }
    
    void Model_DAE::fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot){
        fxdot(t,x->getNVector(),dx->getNVector(),xdot->getNVector());
    }
    
    /** residual function of the DAE
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     */
    void Model_DAE::fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot) {
        computeX_pos(x);
        fw(t,x);
        N_VConst(0.0,xdot);
        fxdot(N_VGetArrayPointer(xdot),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                     N_VGetArrayPointer(dx),w.data());
    }
    
    /** diagonalized Jacobian (for preconditioning)
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @return status flag indicating successful execution
     **/
    void Model_DAE::fJDiag(realtype t, AmiVector *JDiag, realtype cj, AmiVector *x,
                          AmiVector *dx) {
        computeX_pos(x->getNVector());
        fdwdx(t,x_pos.getNVector());
        memset(JDiag->data(),0.0,sizeof(realtype)*nx);
        fJDiag(JDiag->data(),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
               0.0,dx->data(),w.data(),dwdx.data());
        if(!checkFinite(nx,JDiag->data(),"Jacobian"))
            throw AmiException("Evaluation of fJDiag failed!");
    }
    
    /** Sensitivity of dx/dt wrt model parameters p
     * @param t timepoint 
     * @param x Vector with the states 
     * @param dx Vector with the derivative states
     * @return status flag indicating successful execution
     */
     void Model_DAE::fdxdotdp(const realtype t, const N_Vector x, const N_Vector dx) {
         std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
         computeX_pos(x);
         fdwdp(t,x_pos.getNVector());
         for(int ip = 0; ip < nplist(); ip++)
             fdxdotdp(&dxdotdp.at(nx*ip),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                      plist_[ip],N_VGetArrayPointer(dx),w.data(),dwdp.data());
    }
    
    /**
     * @brief Mass matrix for DAE systems
     * @param t timepoint
     * @param x Vector with the states
     */
     void Model_DAE::fM(realtype t, const N_Vector x) {
         std::fill(M.begin(),M.end(),0.0);
         computeX_pos(x);
         fM(M.data(),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data());
    }
    
    std::unique_ptr<Solver> Model_DAE::getSolver() {
        return std::unique_ptr<Solver>(new amici::IDASolver());
    }
    
    /** Jacobian of xBdot with respect to adjoint state xB
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param JB Matrix to which the Jacobian will be written
     **/
    void Model_DAE::fJB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, DlsMat JB) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        SetToZero(JB);
        fJB(JB->data,t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                   cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                   w.data(),dwdx.data());
    }
    
    /** JB in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param JB Matrix to which the Jacobian will be written
     */
    void Model_DAE::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, SlsMat JB) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        SparseSetMatToZero(JB);
        fJSparseB(JB,t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                         cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                         w.data(),dwdx.data());
    }
    
    /** Matrix vector product of JB with a vector v (for iterative solvers)
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be
     *written
     * @param cj scalar in Jacobian (inverse stepsize)
     **/
    void Model_DAE::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                        N_Vector vB, N_Vector JvB, realtype cj) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        N_VConst(0.0,JvB);
        fJvB(N_VGetArrayPointer(JvB),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                    cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                    N_VGetArrayPointer(vB),w.data(),dwdx.data());
    }
    
    /** Right hand side of differential equation for adjoint state xB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     */
    void Model_DAE::fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                          N_Vector dxB, N_Vector xBdot) {
        computeX_pos(x);
        fdwdx(t,x_pos.getNVector());
        N_VConst(0.0,xBdot);
        fxBdot(N_VGetArrayPointer(xBdot),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                      N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                      w.data(),dwdx.data());
    }
    
    /** Right hand side of integral equation for quadrature states qB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void Model_DAE::fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot) {
        computeX_pos(x);
        fdwdp(t,x_pos.getNVector());
        N_VConst(0.0,qBdot);
        realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
        for(int ip = 0; (unsigned)ip < plist_.size(); ip++)
            fqBdot(&qBdot_tmp[ip*nJ],plist_[ip],t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
                          N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                          w.data(),dwdp.data());
    }
    
    void Model_DAE::fsxdot(realtype t, AmiVector *x, AmiVector *dx, int ip,
                           AmiVector *sx, AmiVector *sdx, AmiVector *sxdot) {
        fsxdot(t,x->getNVector(),dx->getNVector(), ip,
               sx->getNVector(),sdx->getNVector(),
               sxdot->getNVector());
    }
    
    /** Right hand side of differential equation for state sensitivities sx
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sdx Vector with the derivative state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     */
    void Model_DAE::fsxdot(realtype t, N_Vector x, N_Vector dx, int ip,
                          N_Vector sx, N_Vector sdx, N_Vector sxdot) {
        computeX_pos(x);
        if(ip == 0) { // we only need to call this for the first parameter index will be the same for all remaining
            fM(t,x_pos.getNVector());
            fdxdotdp(t,x_pos.getNVector(),dx);
            fJSparse(t,0.0,x_pos.getNVector(),dx,J);// also calls dwdx & dx
        }
        N_VConst(0.0,sxdot);
        fsxdot(N_VGetArrayPointer(sxdot),t,x_pos.data(), unscaledParameters.data(),fixedParameters.data(),h.data(),
               plist_[ip],N_VGetArrayPointer(dx),N_VGetArrayPointer(sx),N_VGetArrayPointer(sdx),
               w.data(),dwdx.data(),J->data,M.data(),&dxdotdp.at(ip*nx));
    }
}
