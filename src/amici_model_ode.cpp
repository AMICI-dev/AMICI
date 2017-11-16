#include <include/amici_solver_cvodes.h>
#include <include/amici_model_ode.h>
#include <include/udata.h>

namespace amici {
    
    void Model_ODE::fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                          AmiVector *xdot, DlsMat J) {
        CVodeSolver::fJ(J->N*J->M,t, x->getNVector(), xdot->getNVector(),
           J, this, nullptr, nullptr, nullptr);
        
    }
    
    void Model_ODE::fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                        AmiVector *xdot, SlsMat J){
        CVodeSolver::fJSparse(t,x->getNVector(),xdot->getNVector(),J,this,nullptr,nullptr,nullptr);
    }
    
    void Model_ODE::fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                        AmiVector *v, AmiVector *Jv, realtype cj){
        CVodeSolver::fJv(v->getNVector(),Jv->getNVector(),t,x->getNVector(),xdot->getNVector(),
                    this,nullptr);
    }
    
    void Model_ODE::froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root){
        CVodeSolver::froot(t,x->getNVector(),root,this);
    }
    
    void Model_ODE::fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot){
        CVodeSolver::fxdot(t,x->getNVector(),xdot->getNVector(),this);
    }
    
    /** diagonalized Jacobian (for preconditioning)
     * @param[in] t timepoint
     * @param[out] JDiag Vector to which the Jacobian diagonal will be written
     * @param[in] cj scaling factor, inverse of the step size
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] user_data object with user input @type UserData
     **/
    int Model_ODE::fJDiag(realtype t, AmiVector *JDiag, realtype cj, AmiVector *x,
                          AmiVector *dx) {
        fdwdx(t,x->getNVector());
        memset(JDiag->data(),0.0,sizeof(realtype)*nx);
        if(!model_JDiag(JDiag->data(),t,x->data(),p.data(),k.data(),h.data(),w.data(),dwdx.data()))
            return AMICI_ERROR;
        return CVodeSolver::checkVals(nx,JDiag->data(),"Jacobian");
    }
    
    /** Sensitivity of dx/dt wrt model parameters p
     * @param[in] t timepoint @type realtype
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states
     * @param[in] user_data pointer to temp data object
     */
    int Model_ODE::fdxdotdp(const realtype t, const N_Vector x) {
        std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
        fdwdp(t,x);
        for(int ip = 1; ip < nplist; ip++)
            if(!model_dxdotdp(dxdotdp.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),
                          plist[ip],w.data(),dwdp.data()))
                return AMICI_ERROR;
        return AMICI_SUCCESS;
    }
    
    Solver *Model_ODE::getSolver() {
        return new amici::CVodeSolver();
    }
    
}
