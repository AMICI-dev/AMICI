#ifndef AMICI_MODEL_DAE_H
#define AMICI_MODEL_DAE_H

#include <include/amici_model.h>
#include <include/amici_defines.h>
#include <include/amici_exception.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <vector>

namespace amici {
    
    class UserData;
    class ExpData;
    
    /**
     * @brief The Model class represents an AMICI ODE model.
     * The model does not contain any data, but represents the state
     * of the model at a specific time t. The states must not always be
     * in sync, but may be updated asynchroneously. 
     */
    class Model_DAE : public Model {
    public:
        /** default constructor */
        Model_DAE() : Model() {}
        
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
        Model_DAE(const int np, const int nx, const int nxtrue, const int nk,
              const int ny, const int nytrue, const int nz, const int nztrue,
              const int ne, const int nJ, const int nw, const int ndwdx,
              const int ndwdp, const int nnz, const int ubw, const int lbw,
              const AMICI_o2mode o2mode)
        : Model(np,nx,nxtrue,nk,ny,nytrue,nz,nztrue,ne,nJ,nw,ndwdx,ndwdp,nnz,ubw,lbw,o2mode){}
       

        static int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx,
                      N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1,
                      N_Vector tmp2, N_Vector tmp3);
        
        /** model specific implementation for fJ
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] dx Vector with the derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         **/
        virtual void model_J(realtype *J, const realtype t, const realtype *x, const double *p, const double *k,
                            const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) = 0;

        static int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                       N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                       N_Vector tmp2B, N_Vector tmp3B);
        
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
        virtual int model_JB(realtype *J, const realtype t, const realtype *x, const double *p, const double *k,
                             const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return(AMICI_ERROR);
        }

        static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J,
                            void *user_data, N_Vector tmp1, N_Vector tmp2,
                            N_Vector tmp3);
        
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
        virtual int model_JSparse(realtype *J, const realtype t, const realtype *x, const double *p, const double *k,
                            const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) = 0;
        
        static int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                             SlsMat JB, void *user_data, N_Vector tmp1B,
                             N_Vector tmp2B, N_Vector tmp3B);
        
        /** model specific implementation for fJSparseB
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
        virtual int model_JSparseB(realtype *J, const realtype t, const realtype *x, const double *p, const double *k,
                                   const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB,
                                   const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        }
        
        static int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj,
                          N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
        
        
        static int fJBandB(long int NeqBdot, long int mupper, long int mlower,
                           realtype t, realtype cj, N_Vector x, N_Vector dx,
                           N_Vector xB, N_Vector dxB, N_Vector xBdot,
                           DlsMat JB, void *user_data, N_Vector tmp1B,
                           N_Vector tmp2B, N_Vector tmp3B);

        virtual void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx,
                            void *user_data);
        
        static int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv,
                       realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
        
        /** model specific implementation for fJv
         * @param[out] Jv Matrix vector product of J with a vector v
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] dx Vector with the derivative states
         * @param[in] v Vector with which the Jacobian is multiplied
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_Jv(realtype *Jv, const realtype t, const realtype *x, const double *p, const double *k,
                             const realtype cj, const realtype *dx,const realtype *v,
                             const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        }
        
        static int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                        N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                        N_Vector tmpB1, N_Vector tmpB2);
        
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
        virtual int model_JvB(realtype *JvB, const realtype t, const realtype *x, const double *p, const double *k,
                              const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB,
                              const realtype *vB, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
        }
        
        virtual void frootwrap(realtype t, N_Vector x, N_Vector dx, realtype *root,
                               void *user_data){
            froot(t,x,dx,root,user_data);
        }
        
        static int froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                         void *user_data);
        
        /** model specific implementation for froot
         * @param[out] root values of the trigger function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] dx Vector with the derivative states
         * @return status flag indicating successful execution
         **/
        virtual int model_root(realtype *root, const realtype t, const realtype *x, const double *p, const double *k,
                              const realtype *dx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
        }

        static int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                         void *user_data);
        
        /** model specific implementation for fxdot
         * @param[out] xdot residual function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         * @param[in] dx Vector with the derivative states
         * @return status flag indicating successful execution
         **/
        virtual void model_xdot(realtype *xdot, const realtype t, const realtype *x, const double *p, const double *k,
                                const realtype *dx, const realtype *w) = 0;
        
        static int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                          N_Vector dxB, N_Vector xBdot, void *user_data);
        
        /** model specific implementation for fxBdot
         * @param[out] xBdot adjoint residual function
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] xB Vector with the adjoint states
         * @param[in] dx Vector with the derivative states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_xBdot(realtype *xBdot, const realtype t, const realtype *x, const double *p, const double *k,
                               const realtype *xB, const realtype *dx, const realtype *dxB,
                               const realtype *w, const realtype *dwdx) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
            
        }

        static int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                          void *user_data);
        
        /** model specific implementation for fqBdot
         * @param[out] qBdot adjoint quadrature equation
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] xB Vector with the adjoint states
         * @param[in] dx Vector with the derivative states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         * @return status flag indicating successful execution
         **/
        virtual int model_qBdot(realtype *qBdot, const realtype t, const realtype *x, const double *p, const double *k,
                                const realtype *xB, const realtype *dx, const realtype *dxB,
                                const realtype *w, const realtype *dwdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
            
        }
        
        void fdxdotdp(const realtype t, const N_Vector x, const N_Vector dx, const UserData *udata);
        
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
        virtual int model_dxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                  const int ip, const realtype *dx, const realtype *w, const realtype *dwdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        };

        static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                          N_Vector *sx, N_Vector *sdx, N_Vector *sxdot, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
        
        /** model specific implementation of fsxdot
         * @param[out] sxdot sensitivity rhs
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] ip parameter index
         * @param[in] dx Vector with the derivative states
         * @param[in] sx Vector with the state sensitivities
         * @param[in] sdx Vector with the derivative state sensitivities
         * @param[in] w vector with helper variables
         * @param[in] dwdx derivative of w wrt x
         */
        virtual int model_sxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                const int ip, const realtype *dx, const realtype *sx, const realtype *sdx,
                                const realtype *w, const realtype *dwdx, const realtype *M, const realtype *J,
                                const realtype *dxdotdp) {
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR;
        };
        
        void fM(realtype t, const N_Vector x, const UserData *udata);
        
        /** model specific implementation of fM
         * @param[out] M mass matrix
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         */
        virtual void model_M(realtype *M, const realtype t, const realtype *x, const realtype *p,
                             const realtype *k) {};

    };
    
} // namespace amici

#endif // MODEL_H
