#include <include/amici_solver_idas.h>
#include <include/amici_model_dae.h>

namespace amici {

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
     * @param[in] t timepoint
     * @param[out] JDiag Vector to which the Jacobian diagonal will be written
     * @param[in] cj scaling factor, inverse of the step size
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] user_data object with user input @type UserData
     **/
    int Model_DAE::fJDiag(realtype t, AmiVector *JDiag, realtype cj, AmiVector *x,
                          AmiVector *dx) {
        fdwdx(t,x->getNVector());
        memset(JDiag->data(),0.0,sizeof(realtype)*nx);
        if(!model_JDiag(JDiag->data(),t,x->data(),p.data(),k.data(),
                               dx->data(),w.data(),dwdx.data()))
            return AMICI_ERROR;
        return IDASolver::checkVals(nx,JDiag->data(),"Jacobian");
    }
    
    /** Sensitivity of dx/dt wrt model parameters p
     * @param[in] t timepoint @type realtype
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states
     * @param[in] user_data pointer to temp data object
     */
     int Model_DAE::fdxdotdp(const realtype t, const N_Vector x, const N_Vector dx) {
        std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
        fdwdp(t,x);
        for(int ip = 1; ip < nplist; ip++)
            if(!model_dxdotdp(dxdotdp.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),
                          plist[ip],N_VGetArrayPointer(dx),w.data(),dwdp.data()))
                return AMICI_ERROR;
        return AMICI_SUCCESS;
    }
    
    /**
     * @brief Mass matrix for DAE systems
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] udata object with model specifications
     */
     void Model_DAE::fM(realtype t, const N_Vector x) {
        std::fill(M.begin(),M.end(),0.0);
        model_M(M.data(),t,N_VGetArrayPointer(x),p.data(),k.data());
    }
    
    Solver *Model_DAE::getSolver() {
        return new amici::IDASolver();
    }
}
