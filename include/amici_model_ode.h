#ifndef AMICI_MODEL_ODE_H
#define AMICI_MODEL_ODE_H

#include <include/amici_model.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <vector>

namespace amici {
    extern msgIdAndTxtFp warnMsgIdAndTxt;
    
    class CVodeSolver;
    
    /**
     * @brief The Model class represents an AMICI ODE model.
     * The model does not contain any data, but represents the state
     * of the model at a specific time t. The states must not always be
     * in sync, but may be updated asynchroneously. 
     */
    class Model_ODE : public Model {
    public:
        /** default constructor */
        Model_ODE() : Model() {}
        
        /** constructor with model dimensions
         * @param np number of parameters
         * @param nx number of state variables
         * @param nxtrue number of state variables of the non-augmented model
         * @param nk number of constants
         * @param ny number of observables
         * @param nytrue number of observables of the non-augmented model
         * @param nz number of event observables
         * @param nztrue number of event observables of the non-augmented model
         * @param ne number of events
         * @param nJ number of objective functions
         * @param nw number of repeating elements
         * @param ndwdx number of nonzero elements in the x derivative of the
         * repeating elements
         * @param ndwdp number of nonzero elements in the p derivative of the
         * repeating elements
         * @param nnz number of nonzero elements in Jacobian
         * @param ubw upper matrix bandwidth in the Jacobian
         * @param lbw lower matrix bandwidth in the Jacobian
         * @param o2mode second order sensitivity mode
         */
        Model_ODE(const int np, const int nx, const int nxtrue, const int nk,
                  const int ny, const int nytrue, const int nz, const int nztrue,
                  const int ne, const int nJ, const int nw, const int ndwdx,
                  const int ndwdp, const int nnz, const int ubw, const int lbw,
                  const AMICI_o2mode o2mode, const std::vector<realtype> p,
                  const std::vector<realtype> k, const std::vector<int> plist,
                  const std::vector<realtype> idlist, const std::vector<int> z2event)
        : Model(np,nx,nxtrue,nk,ny,nytrue,nz,nztrue,ne,nJ,nw,ndwdx,ndwdp,nnz,ubw,lbw,o2mode,p,k,plist,idlist,z2event){}
        
        virtual void fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                        AmiVector *xdot, DlsMat J) override;
        
        /** model specific implementation for fJ
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         **/
        virtual void model_J(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                             const realtype *w, const realtype *dwdx) = 0;
        
        /** model specific implementation for fJB
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] xB Vector with the adjoint states
         * @param[in] dx Vector with the derivative states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_JB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                             const realtype *xB, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return(AMICI_ERROR);
        }
        
        virtual void fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                                  AmiVector *xdot, SlsMat J) override;
        
        /** model specific implementation for fJSparse
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] dx Vector with the derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_JSparse(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                  const realtype *w, const realtype *dwdx) = 0;
        
        /** model specific implementation for fJSparseB
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] xB Vector with the adjoint states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_JSparseB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                   const realtype *xB, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        }
        
        /** model specific implementation for fJSparseB
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_JDiag(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        }
        
        virtual int fJDiag(realtype t, AmiVector *Jdiag, realtype cj, AmiVector *x,
                                AmiVector *dx) override;
        
        virtual void fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                             AmiVector *v, AmiVector *nJv, realtype cj) override;
        
        /** model specific implementation for fJv
         * @param[out] Jv Matrix vector product of J with a vector v
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] v Vector with which the Jacobian is multiplied
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_Jv(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                             const realtype *v, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        }
        
        /** model specific implementation for fJvB
         * @param[out] JvB Matrix vector product of JB with a vector v
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] xB Vector with the adjoint states
         * @param[in] dx Vector with the derivative states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] vB Vector with which the Jacobian is multiplied
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_JvB(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                              const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
        }
        
        virtual void froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root) override;
        
        /** model specific implementation for froot
         * @param[out] root values of the trigger function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @return status flag indicating successful execution
         **/
        virtual int model_root(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
        }
        
        virtual void fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot) override;
        
        /** model specific implementation for fxdot
         * @param[out] xdot residual function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         * @return status flag indicating successful execution
         **/
        virtual void model_xdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                const realtype *w) = 0;
        
        /** model specific implementation for fxBdot
         * @param[out] xBdot adjoint residual function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] xB Vector with the adjoint states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_xBdot(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                const realtype *xB, const realtype *w, const realtype *dwdx) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
            
        }
        
        /** model specific implementation for fqBdot
         * @param[out] qBdot adjoint quadrature equation
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] xB Vector with the adjoint states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_qBdot(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                const realtype *xB, const realtype *w, const realtype *dwdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
            
        }
        
        int fdxdotdp(const realtype t, const N_Vector x);
        
        virtual int fdxdotdp(realtype t, AmiVector *x, AmiVector *dx) override {
            return fdxdotdp(t,x->getNVector());
        };
        
        /** model specific implementation of fdxdotdp
         * @param[out] dxdotdp partial derivative xdot wrt p
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] ip parameter index
         * @param[in] dx Vector with the derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdp derivative of w wrt p
         */
        virtual int model_dxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                  const int ip, const realtype *w, const realtype *dwdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        };
        
        /** model specific implementation of fsxdot
         * @param[out] sxdot sensitivity rhs
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] ip parameter index
         * @param[in] sx Vector with the state sensitivities
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         */
        virtual int model_sxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J,
                                const realtype *dxdotdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        };
        
        virtual Solver *getSolver() override;
        
        friend class CVodeSolver;
    };
    
} // namespace amici

#endif // MODEL_H
