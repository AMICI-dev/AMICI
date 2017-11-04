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
         * @param[in] N number of state variables @type long_int
         * @param[in] t timepoint @type realtype
         * @param[in] cj scaling factor, inverse of the step size (only DAE) @type
         *realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xdot Vector with the right hand side @type N_Vector
         * @param[out] J Matrix to which the Jacobian will be written @type DlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1 temporary storage vector @type N_Vector
         * @param[in] tmp2 temporary storage vector @type N_Vector
         * @param[in] tmp3 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         **/
        virtual int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3) {
            UserData udata = (UserData) user_data;
            return model_J(J,t,cj,x,p,k,dx,xdot)
        }

        /** model specific implementation of fJ
        virtual int model_J(...) = 0;

                /** Jacobian of xBdot with respect to adjoint state xB
         * @param[in] NeqBdot number of adjoint state variables @type long_int
         * @param[in] t timepoint @type realtype
         * @param[in] cj scaling factor, inverse of the step size (only DAE) @type
         *realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[in] xBdot Vector with the adjoint right hand side @type N_Vector
         * @param[out] JB Matrix to which the Jacobian will be written @type DlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1B temporary storage vector @type N_Vector
         * @param[in] tmp2B temporary storage vector @type N_Vector
         * @param[in] tmp3B temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         **/
        virtual int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                         N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
            throw AmiException("Missing function implementation!");
        }

        /** diagonalized Jacobian (for preconditioning)
         * @param[in] t timepoint @type realtype
         * @param[out] JDiag Vector to which the Jacobian diagonal will be written
         *@type NVector
         * @param[in] cj scaling factor, inverse of the step size (only DAE) @type
         *realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] user_data object with model specifications @type TempData
         **/
        virtual void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx,
                            void *user_data) {
            throw AmiException("Missing function implementation!");
        }

        /** Matrix vector product of J with a vector v (for iterative solvers)
         * @param[in] t timepoint @type realtype
         * @param[in] cj scaling factor, inverse of the step size (only DAE) @type
         *realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xdot Vector with the right hand side @type N_Vector
         * @param[in] v Vector with which the Jacobian is multiplied @type N_Vector
         * @param[out] Jv Vector to which the Jacobian vector product will be
         *written @type N_Vector
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1 temporary storage vector @type N_Vector
         * @param[in] tmp2 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         **/
        virtual int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv,
                         realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
            throw AmiException("Missing function implementation!");
        }

        /** Event trigger function for events
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         * N_Vector
         * @param[out] root array with root function values @type realtype
         * @param[in] user_data object with model specifications @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual int froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                          void *user_data) {
            throw AmiException("Missing function implementation!");
        }

                /** Right hand side of differential equation for states x
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         * N_Vector
         * @param[out] xdot Vector with the right hand side @type N_Vector
         * @param[in] user_data pointer to temp data object @type TempDat
         * @return status flag indicating successful execution @type int
         */
        virtual int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                          void *user_data) {
            throw AmiException("Missing function implementation!");
        }
        
        /** Right hand side of differential equation for adjoint state xB
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         * N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[out] xBdot Vector with the adjoint right hand side @type N_Vector
         * @param[in] user_data pointer to temp data object @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                           N_Vector dxB, N_Vector xBdot, void *user_data) {
            throw AmiException("Missing function implementation!");
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
            throw AmiException("Missing function implementation!");
        }
        
        /** Sensitivity of dx/dt w.r.t. model parameters p
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         * N_Vector
         * @param[in] user_data pointer to temp data object @type TempData
         */
        virtual void fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data) {
            throw AmiException("Missing function implementation!");
        }

        /** Right hand side of differential equation for state sensitivities sx
         * @param[in] Ns number of parameters @type int
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xdot Vector with the right hand side @type N_Vector
         * @param[in] ip parameter index @type int
         * @param[in] sx Vector with the state sensitivities @type N_Vector
         * @param[in] sdx Vector with the derivative state sensitivities (only DAE)
         * @type N_Vector
         * @param[out] sxdot Vector with the sensitivity right hand side @type
         * N_Vector
         * @param[in] user_data pointer to temp data object @type TempData
         * @param[in] tmp1 temporary storage vector @type N_Vector
         * @param[in] tmp2 temporary storage vector @type N_Vector
         * @param[in] tmp3 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         */
        virtual int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot, int ip,
                            N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
            throw AmiException("Missing function implementation!");
        }
        
        /** J in sparse form (for sparse solvers from the SuiteSparse Package)
         * @param[in] t timepoint @type realtype
         * @param[in] cj scalar in Jacobian (inverse stepsize, only DAE) @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xdot Vector with the right hand side @type N_Vector
         * @param[out] J Matrix to which the Jacobian will be written @type SlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1 temporary storage vector @type N_Vector
         * @param[in] tmp2 temporary storage vector @type N_Vector
         * @param[in] tmp3 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         */
        virtual int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J,
                              void *user_data, N_Vector tmp1, N_Vector tmp2,
                              N_Vector tmp3) {
            throw AmiException("Missing function implementation!");
        }
        
        /** J in banded form (for banded solvers)
         * @param[in] N number of states @type long int
         * @param[in] mupper upper matrix bandwidth @type long int
         * @param[in] mlower lower matrix bandwidth @type long int
         * @param[in] t timepoint @type realtype
         * @param[in] cj scalar in Jacobian (inverse stepsize, only DAE) @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xdot Vector with the right hand side @type N_Vector
         * @param[out] J Matrix to which the Jacobian will be written @type DlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1 temporary storage vector @type N_Vector
         * @param[in] tmp2 temporary storage vector @type N_Vector
         * @param[in] tmp3 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         */
        virtual int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj,
                            N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
            throw AmiException("Missing function implementation!");
        }
        
        /** JB in banded form (for banded solvers)
         * @param[in] NeqBdot number of states @type long int
         * @param[in] mupper upper matrix bandwidth @type long int
         * @param[in] mlower lower matrix bandwidth @type long int
         * @param[in] t timepoint @type realtype
         * @param[in] cj scalar in Jacobian (inverse stepsize, only DAE) @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[in] xBdot Vector with the adjoint right hand side @type N_Vector
         * @param[out] JB Matrix to which the Jacobian will be written @type DlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1B temporary storage vector @type N_Vector
         * @param[in] tmp2B temporary storage vector @type N_Vector
         * @param[in] tmp3B temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         */
        virtual int fJBandB(long int NeqBdot, long int mupper, long int mlower,
                             realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                             DlsMat JB, void *user_data, N_Vector tmp1B,
                             N_Vector tmp2B, N_Vector tmp3B) {
            throw AmiException("Missing function implementation!");
        }
        
        /** Matrix vector product of JB with a vector v (for iterative solvers)
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[in] xBdot Vector with the adjoint right hand side @type N_Vector
         * @param[in] vB Vector with which the Jacobian is multiplied @type
         *N_Vector
         * @param[out] JvB Vector to which the Jacobian vector product will be
         *written @type N_Vector
         * @param[in] cj scalar in Jacobian (inverse stepsize, only DAE) @type realtype
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmpB1 temporary storage vector @type N_Vector
         * @param[in] tmpB2 temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         **/
        virtual int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                          N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                          N_Vector tmpB1, N_Vector tmpB2) {
            throw AmiException("Missing function implementation!");
        }
        
        /** JB in sparse form (for sparse solvers from the SuiteSparse Package)
         * @param[in] t timepoint @type realtype
         * @param[in] cj scalar in Jacobian (inverse stepsize, only DAE) @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] xB Vector with the adjoint states @type N_Vector
         * @param[in] dxB Vector with the adjoint derivative states (only DAE)
         * @type N_Vector
         * @param[in] xBdot Vector with the adjoint right hand side @type N_Vector
         * @param[out] JB Matrix to which the Jacobian will be written @type DlsMat
         * @param[in] user_data object with model specifications @type TempData
         * @param[in] tmp1B temporary storage vector @type N_Vector
         * @param[in] tmp2B temporary storage vector @type N_Vector
         * @param[in] tmp3B temporary storage vector @type N_Vector
         * @return status flag indicating successful execution @type int
         */
        virtual int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                               SlsMat JB, void *user_data, N_Vector tmp1B,
                               N_Vector tmp2B, N_Vector tmp3B) {
            throw AmiException("Missing function implementation!");
        }

                /**
         * @brief Mass matrix for DAE systems (only DAE)
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] user_data object with model specifications @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual void fM(realtype t, N_Vector x, N_Vector dx, void *user_data) {
            throw AmiException("Missing function implementation!"); 
        }
        
        /**
         * @brief jacobian of the right hand side (only DAE)
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] user_data object with model specifications @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual void fdfdx(realtype t, N_Vector x, N_Vector dx, void *user_data) {
            throw AmiException("Missing function implementation!"); 
        }
        
        /**
         * @brief Recurring terms in xdot
         * @param[in] udata object with user input
         */
        virtual void fw(const UserData *udata) {
            throw AmiException("Missing function implementation!");
        }
        
        /**
         * @brief Recurring terms in xdot, parameter derivative
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] user_data object with model specifications @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual void fdwdp(realtype t, N_Vector x, N_Vector dx, void *user_data) {
            throw AmiException("Missing function implementation!");
        }
        
        /**
         * @brief Recurring terms in xdot, state derivative
         * @param[in] t timepoint @type realtype
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] dx Vector with the derivative states (only DAE) @type
         *N_Vector
         * @param[in] user_data object with model specifications @type TempData
         * @return status flag indicating successful execution @type int
         */
        virtual void fdwdx(realtype t, N_Vector x, N_Vector dx, void *user_data) {
            throw AmiException("Missing function implementation!");
        }

    };
    
} // namespace amici

#endif // MODEL_H
