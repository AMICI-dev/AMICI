#include <include/amici_solver_cvodes.h>
#include <include/amici_model_ode.h>
#include <include/udata.h>

namespace amici {
    
    void Model_ODE::fJ(realtype t, realtype cj, AmiVector x, AmiVector dx,
                          AmiVector xdot, DlsMat J, CVodeSolver *solver) {
        solver->fJ(J->N*J->M,t, x.getNVector(), xdot.getNVector(),
           J, this, nullptr, nullptr, nullptr);
        
    }
    
    void Model_ODE::fJSparse(realtype t, realtype cj, AmiVector x, AmiVector dx,
                        AmiVector xdot, SlsMat J, CVodeSolver *solver){
        solver->fJSparse(t,x.getNVector(),xdot.getNVector(),J,this,nullptr,nullptr,nullptr);
    }
    
    void Model_ODE::froot(realtype t, AmiVector x, AmiVector dx, realtype *root,
                          CVodeSolver *solver){
        solver->froot(t,x.getNVector(),root,this);
    }
    
    void Model_ODE::fxdot(realtype t, AmiVector x, AmiVector dx, AmiVector xdot,
                          CVodeSolver *solver){
        solver->fxdot(t,x.getNVector(),xdot.getNVector(),this);
    }
    
    /** Sensitivity of dx/dt wrt model parameters p
     * @param[in] t timepoint @type realtype
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states
     * @param[in] user_data pointer to temp data object
     */
    void Model_ODE::fdxdotdp(const realtype t, const N_Vector x) {
        std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
        fdwdp(t,x);
        for(int ip = 1; ip < nplist; ip++)
            model_dxdotdp(dxdotdp.data(),t,N_VGetArrayPointer(x),p.data(),k.data(),h.data(),
                          plist[ip],w.data(),dwdp.data());
    }
    
}
