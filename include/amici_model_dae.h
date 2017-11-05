#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

#include <include/amici.h>
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
    class Model_DAE {
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
        : Model(np,nx,nxtru,nk,ny,nytrue,nz,nztrue,ne,nJ,nw,ndwdx,ndwdp,nnz,ubw,lbw,o2mode){}
       
        /** Jacobian of xdot with respect to states x
         * @param[in] N number of state variables
         * @param[in] t timepoint
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xdot Vector with the right hand side
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1 temporary storage vector
         * @param[in] tmp2 temporary storage vector
         * @param[in] tmp3 temporary storage vector
         * @return status flag indicating successful execution
         **/
        virtual int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(J->data,0.0,sizeof(realtype)*N)
            model_J(J->data,t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                    cj,N_VGetArrayPointer(dx),w.data(),dwdx.data());
            return checkVals(N,J->data,"Jacobian");
        }
        
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

        /** Jacobian of xBdot with respect to adjoint state xB
         * @param[in] NeqBdot number of adjoint state variables
         * @param[in] t timepoint
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xB Vector with the adjoint states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] xBdot Vector with the adjoint right hand side
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1B temporary storage vector
         * @param[in] tmp2B temporary storage vector
         * @param[in] tmp3B temporary storage vector
         * @return status flag indicating successful execution
         **/
        virtual int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                         N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(JB->data,0.0,sizeof(realtype)*NeqBdot)
            if(!model_JB(JB->data,t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                         cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                         w.data(),dwdx.data()))
                return AMICI_ERROR;
            return checkVals(N,J->data,"Jacobian");
        }
        
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
        
        /** J in sparse form (for sparse solvers from the SuiteSparse Package)
         * @param[in] t timepoint
         * @param[in] cj scalar in Jacobian (inverse stepsize)
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         *N_Vector
         * @param[in] xdot Vector with the right hand side
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1 temporary storage vector
         * @param[in] tmp2 temporary storage vector
         * @param[in] tmp3 temporary storage vector
         * @return status flag indicating successful execution
         */
        virtual int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J,
                             void *user_data, N_Vector tmp1, N_Vector tmp2,
                             N_Vector tmp3) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(J->data,0.0,sizeof(realtype)*J->nnz)
            model_JSparse(J->data,t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                          cj,N_VGetArrayPointer(dx),w.data(),dwdx.data());
            return checkVals(J->nnz,J->data,"Jacobian");
        }
        
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
        
        
        /** JB in sparse form (for sparse solvers from the SuiteSparse Package)
         * @param[in] t timepoint
         * @param[in] cj scalar in Jacobian
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xB Vector with the adjoint states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] xBdot Vector with the adjoint right hand side
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1B temporary storage vector
         * @param[in] tmp2B temporary storage vector
         * @param[in] tmp3B temporary storage vector
         * @return status flag indicating successful execution
         */
        virtual int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                              SlsMat JB, void *user_data, N_Vector tmp1B,
                              N_Vector tmp2B, N_Vector tmp3B) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(JB->data,0.0,sizeof(realtype)*JB->nnz)
            if(!model_JSparseB(JB->data,t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                               cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                               w.data(),dwdx.data()))
                return AMICI_ERROR;
            return checkVals(N,J->data,"Jacobian");
        }
        
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
        
        /** J in banded form (for banded solvers)
         * @param[in] N number of states
         * @param[in] mupper upper matrix bandwidth
         * @param[in] mlower lower matrix bandwidth
         * @param[in] t timepoint
         * @param[in] cj scalar in Jacobian (inverse stepsize)
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xdot Vector with the right hand side
         * @param[out] J Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1 temporary storage vector
         * @param[in] tmp2 temporary storage vector
         * @param[in] tmp3 temporary storage vector
         * @return status flag indicating successful execution
         */
        virtual int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj,
                           N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
            return fJ(N,t,cj,x,dx,xdot,J,user_data,tmp1,tmp2,tmp3);
        }
        
        
        
        /** JB in banded form (for banded solvers)
         * @param[in] NeqBdot number of states
         * @param[in] mupper upper matrix bandwidth
         * @param[in] mlower lower matrix bandwidth
         * @param[in] t timepoint
         * @param[in] cj scalar in Jacobian (inverse stepsize)
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xB Vector with the adjoint states
         * @param[in] dxB Vector with the adjoint derivative states
         * @param[in] xBdot Vector with the adjoint right hand side
         * @param[out] JB Matrix to which the Jacobian will be written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1B temporary storage vector
         * @param[in] tmp2B temporary storage vector
         * @param[in] tmp3B temporary storage vector
         * @return status flag indicating successful execution
         */
        virtual int fJBandB(long int NeqBdot, long int mupper, long int mlower,
                            realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                            DlsMat JB, void *user_data, N_Vector tmp1B,
                            N_Vector tmp2B, N_Vector tmp3B) {
            return fJB(NeqBdot,t,cj,x,dx,xB,dxB,xBdot,JB,user_data,tmp1B,tmp2B,tmp3B);
        }


        /** diagonalized Jacobian (for preconditioning)
         * @param[in] t timepoint
         * @param[out] JDiag Vector to which the Jacobian diagonal will be written
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] user_data object with user input @type UserData
         **/
        virtual void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx,
                            void *user_data) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(N_VGetArrayPointer(JDiag),0.0,sizeof(realtype)*nx)
            if(!model_JDiag(N_VGetArrayPointer(JDiag),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                               cj,N_VGetArrayPointer(dx),w.data(),dwdx.data()))
                return AMICI_ERROR;
            return checkVals(nx,N_VGetArrayPointer(JDiag),"Jacobian");
        }

        /** Matrix vector product of J with a vector v (for iterative solvers)
         * @param[in] t timepoint @type realtype
         * @param[in] cj scaling factor, inverse of the step size
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xdot Vector with the right hand side
         * @param[in] v Vector with which the Jacobian is multiplied
         * @param[out] Jv Vector to which the Jacobian vector product will be
         *written
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1 temporary storage vector
         * @param[in] tmp2 temporary storage vector
         * @return status flag indicating successful execution
         **/
        virtual int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv,
                         realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(N_VGetArrayPointer(Jv),0.0,sizeof(realtype)*nx)
            if(!model_Jv(N_VGetArrayPointer(Jv),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                            cj,N_VGetArrayPointer(dx),N_VGetArrayPointer(v),w.data(),dwdx.data()))
                return AMICI_ERROR;
            return checkVals(nx,N_VGetArrayPointer(Jv),"Jacobian");
        }
        
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
        
        /** Matrix vector product of JB with a vector v (for iterative solvers)
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xB Vector with the adjoint states
         * @param[in] dxB Vector with the adjoint derivative states
         * @type N_Vector
         * @param[in] xBdot Vector with the adjoint right hand side
         * @param[in] vB Vector with which the Jacobian is multiplied 
         * @param[out] JvB Vector to which the Jacobian vector product will be
         *written
         * @param[in] cj scalar in Jacobian (inverse stepsize)
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmpB1 temporary storage vector
         * @param[in] tmpB2 temporary storage vector
         * @return status flag indicating successful execution
         **/
        virtual int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                         N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                         N_Vector tmpB1, N_Vector tmpB2) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(N_VGetArrayPointer(JvB),0.0,sizeof(realtype)*nx)
            if(!model_JvB(N_VGetArrayPointer(JvB),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                          cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                          N_VGetArrayPointer(vB),w.data(),dwdx.data()))
                return AMICI_ERROR;
            return checkVals(nx,N_VGetArrayPointer(JvB),"Jacobian");
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
        virtual int model_JvB(realtype *JvB, const realtype t, const realtype *x, const double *p, const double *k,
                              const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB,
                              const realtype *vB, const realtype *w, const realtype *dwdx){
            warnMsgIdAndTxt("AMICI:mex","Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
            return AMICI_ERROR; // not implemented
        }
        

        /** Event trigger function for events
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[out] root array with root function values
         * @param[in] user_data object with user input @type UserData
         * @return status flag indicating successful execution
         */
        virtual int froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                          void *user_data) {
            UserData udata = (UserData) user_data;
            memset(root,0.0,sizeof(realtype)*ne)
            if(!model_root(root,t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                           N_VGetArrayPointer(dx))) {
                return AMICI_ERROR;
            }
            return checkVals(ne,root,"root function");
            
        }
        
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

        /** residual function of the DAE
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states 
         * @param[out] xdot Vector with the right hand side
         * @param[in] user_data object with user input @type UserData
         * @return status flag indicating successful execution
         */
        virtual int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                          void *user_data) {
            UserData udata = (UserData) user_data;
            fw(t,x,udata);
            memset(N_VGetArrayPointer(xdot),0.0,sizeof(realtype)*nx)
            model_xdot(N_VGetArrayPointer(xdot),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                       N_VGetArrayPointer(dx),w.data());
            return checkVals(nx,xdot,"residual function");
        }
        
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
        
        /** Right hand side of differential equation for adjoint state xB
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states 
         * @param[in] xB Vector with the adjoint states
         * @param[in] dxB Vector with the adjoint derivative states 
         * @param[out] xBdot Vector with the adjoint right hand side
         * @param[in] user_data object with user input @type UserData
         * @return status flag indicating successful execution
         */
        virtual int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                           N_Vector dxB, N_Vector xBdot, void *user_data) {
            UserData udata = (UserData) user_data;
            fdwdx(t,x,udata);
            memset(N_VGetArrayPointer(xBdot),0.0,sizeof(realtype)*nx)
            model_xBdot(N_VGetArrayPointer(xBdot),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                        N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                        w.data(),dwdx.data());
            return checkVals(nx,N_VGetArrayPointer(xBdot),"adjoint residual function");
        }
        
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

        /** Right hand side of integral equation for quadrature states qB
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[out] qBdot Vector with the adjoint quadrature right hand side
         * @type N_Vector
         * @param[in] user_data pointer to temp data object @type TempDat
         * @return status flag indicating successful execution @type int
         */
        virtual int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                           void *user_data) {
            UserData udata = (UserData) user_data;
            fdwdp(t,x,udata);
            memset(N_VGetArrayPointer(qBdot),0.0,sizeof(realtype)*nJ*udata->nplist)
            model_qBdot(N_VGetArrayPointer(xBdot),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                        N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                        w.data(),dwdp.data());
            return checkVals(nx,N_VGetArrayPointer(qBdot),"adjoint quadrature function");
        }
        
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
        
        /** Sensitivity of dx/dt wrt model parameters p
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states
         * @param[in] user_data pointer to temp data object
         */
        virtual void fdxdotdp(const realtype t, const N_Vector x, const N_Vector dx, const UserData *udata) {
            UserData udata = (UserData) user_data;
            std::fill(dxdotdp.begin(),dxdotdp.end(),0.0);
            fdwdp(t,x,udata);
            for(int ip = 1; ip < nplist; ip++)
                model_dxdotdp(dxdotdp.data(),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                              udata->plist[ip],N_VGetArrayPointer(dx),w.data(),dwdP.data());
        }
        
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

        /** Right hand side of differential equation for state sensitivities sx
         * @param[in] Ns number of parameters
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states
         * @param[in] xdot Vector with the right hand side
         * @param[in] ip parameter index
         * @param[in] sx Vector with the state sensitivities
         * @param[in] sdx Vector with the derivative state sensitivities
         * @param[out] sxdot Vector with the sensitivity right hand side
         * @param[in] user_data object with user input @type UserData
         * @param[in] tmp1 temporary storage vector
         * @param[in] tmp2 temporary storage vector
         * @param[in] tmp3 temporary storage vector
         * @return status flag indicating successful execution
         */
        virtual int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot, int ip,
                            N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
            UserData udata = (UserData) user_data;
            if(ip == 0) { // we only need to call this for the first parameter index will be the same for all remaining
                fM(t,x,udata);
                if(!fdxdotdp(t,x,dx,udata))
                    return AMICI_ERROR;
                if(!fJ(nx*nx,t,0.0,x,dx,nullptr,J,udata,tmp1,tmp2,tmp3))// also calls dwdx & dx
                    return AMICI_ERROR;
            }
            memset(N_VGetArrayPointer(sxdot),0.0,sizeof(realtype)*nx)
            if(!model_sxdot(N_VGetArrayPointer(sxdot),t,N_VGetArrayPointer(x),udata->p(),udata->k(),
                            udata->plist[ip],N_VGetArrayPointer(dx),N_VGetArrayPointer(sx),N_VGetArrayPointer(sdx),
                            w.data(),dwdx.data(),M.data(),J.data(),dxdotdp.data())) {
                return AMICI_ERROR;
            }
            return checkVals(nx,N_VGetArrayPointer(sxdot),"sensitivity rhs");
        }
        
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
        
        /**
         * @brief Mass matrix for DAE systems
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] dx Vector with the derivative states 
         * @param[in] udata object with model specifications
         */
        virtual void fM(realtype t, const N_Vector x, const UserData *udata) {
            std::fill(M.begin(),M.end(),0.0)
            M_model(M.data(),t,N_VGetArrayPointer(x),udata->p(),udata->k());
        }
        
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
