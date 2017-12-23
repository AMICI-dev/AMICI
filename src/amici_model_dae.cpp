#include <include/amici_solver_idas.h>
#include <include/amici_model_dae.h>

namespace amici {
    
    void Model_DAE::fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                       AmiVector *xdot, DlsMat J) {
        IDASolver::fJ(J->N*J->M,t,cj, x->getNVector(), dx->getNVector(), xdot->getNVector(),
                        J, this, nullptr, nullptr, nullptr);
    }

     void Model_DAE::fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                              AmiVector *xdot, SlsMat J){
        IDASolver::fJSparse(t,cj,x->getNVector(),dx->getNVector(),xdot->getNVector(),J,this,nullptr,nullptr,nullptr);
    }
    
    void Model_DAE::fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                        AmiVector *v, AmiVector *Jv, realtype cj){
        IDASolver::fJv(t,x->getNVector(),dx->getNVector(),xdot->getNVector(),v->getNVector(),
                    Jv->getNVector(),cj,this,nullptr,nullptr);
    }
    
     void Model_DAE::froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root){
        IDASolver::froot(t,x->getNVector(),dx->getNVector(),root,this);
    }
    
     void Model_DAE::fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot){
        IDASolver::fxdot(t,x->getNVector(),dx->getNVector(),xdot->getNVector(),this);
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
        fdwdx(t,x->getNVector());
        memset(JDiag->data(),0.0,sizeof(realtype)*nx);
        fJDiag(JDiag->data(),t,x->data(),p.data(),k.data(),h.data(),
               0.0,dx->data(),w.data(),dwdx.data());
        if(!isFinite(nx,JDiag->data(),"Jacobian"))
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
        fdwdp(t,x);
        for(int ip = 0; ip < nplist(); ip++)
            fdxdotdp(&dxdotdp.at(nx*ip),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),
                     plist[ip],N_VGetArrayPointer(dx),w.data(),dwdp.data());
    }
    
    /**
     * @brief Mass matrix for DAE systems
     * @param t timepoint
     * @param x Vector with the states
     */
     void Model_DAE::fM(realtype t, const N_Vector x) {
        std::fill(M.begin(),M.end(),0.0);
        fM(M.data(),t,N_VGetArrayPointer(x),p.data(),k.data());
    }
    
    std::unique_ptr<Solver> Model_DAE::getSolver() {
        return std::unique_ptr<Solver>(new amici::IDASolver());
    }
}
