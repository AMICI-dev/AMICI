#ifndef MODEL_H
#define MODEL_H

#include <include/amici.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

class UserData;
class ExpData;

/**
 * @brief The Model class represents an AMICI ODE model.
 * The model does not contain any data, its state should never change.
 */
class Model {
  public:
    /** default constructor */
    Model()
        : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
          ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
          o2mode(AMICI_O2MODE_NONE) {}

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
    Model(const int np, const int nx, const int nxtrue, const int nk,
          const int ny, const int nytrue, const int nz, const int nztrue,
          const int ne, const int nJ, const int nw, const int ndwdx,
          const int ndwdp, const int nnz, const int ubw, const int lbw,
          const AMICI_o2mode o2mode)
        : np(np), nk(nk), nx(nx), nxtrue(nxtrue), ny(ny), nytrue(nytrue),
          nz(nz), nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
          nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode) {}

    /**
     * @brief Returns a UserData instance with preset model dimensions
     * @return The UserData instance
     */
    UserData getUserData() const;

    /**
     * @brief Returns a new UserData instance with preset model dimensions
     * @return The UserData instance
     */
    UserData *getNewUserData() const;

    /** Retrieves the solver object
     * @return Solver solver object @type Solver
     */
    virtual Solver *getSolver() { return NULL; }

    /** Initial states
      * @param[out] x0 Vector to which the initial states will be written @type
      *N_Vector
      * @param[in] user_data object with model specifications @type TempData
      * @return status flag indicating successful execution @type int
     **/
    virtual int fx0(N_Vector x0, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Initial value for time derivative of states (only necessary for DAEs)
      * @param[in] x0 Vector with the initial states @type N_Vector
      * @param[out] dx0 Vector to which the initial derivative states will be
      *written (only DAE) @type N_Vector
      * @param[in] user_data object with model specifications @type TempData
      * @return status flag indicating successful execution @type int
     **/
    virtual int fdx0(N_Vector x0, N_Vector dx0, void *user_data) {
        return AMICI_SUCCESS;
    }

    //    virtual int fdx0(const realtype *k, realtype *x) { return
    //    AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Initial value for initial state sensitivities
      * @param[out] sx0 Vector to whith the initial state sensitivities @type
      *N_Vector
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] dx Vector with the derivative states (only DAE) @type
      *N_Vector
      * @param[in] user_data object with model specifications @type TempData
      * @return status flag indicating successful execution @type int
     **/
    virtual int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of derivative initial states sensitivities sdx0 (only
      *necessary for DAEs)
      * @param[out] sdx0 Vector to whith the derivative initial state
      *sensitivities @type N_Vector
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] dx Vector with the derivative states (only DAE) @type
      *N_Vector
      * @param[in] user_data object with model specifications @type TempData
      * @return status flag indicating successful execution @type int
     **/
    virtual int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx,
                      void *user_data) {
        return AMICI_SUCCESS;
    }
    
    /** Initial value for initial state sensitivities
     * @param[out] sx0 Vector to whith the initial state sensitivities @type
     *N_Vector
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states (only DAE) @type
     *N_Vector
     * @param[in] user_data object with model specifications @type TempData
     * @return status flag indicating successful execution @type int
     **/
    virtual int fs2x0(realtype *s2x0, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Sensitivity of derivative initial states sensitivities sdx0 (only
     *necessary for DAEs)
     * @param[out] sdx0 Vector to whith the derivative initial state
     *sensitivities @type N_Vector
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states (only DAE) @type
     *N_Vector
     * @param[in] user_data object with model specifications @type TempData
     * @return status flag indicating successful execution @type int
     **/
    virtual int fs2dx0(realtype *s2dx0, N_Vector x, N_Vector dx,
                      void *user_data) {
        return AMICI_SUCCESS;
    }

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
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
      * @return status flag indicating successful execution @type int
     **/
    virtual int fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx,
                       void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Event root function of events (equal to froot but does not include
     * non-output events)
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int frz(realtype t, int ie, N_Vector x, TempData *tdata,
                    ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of rz, total derivative
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] sx Vector with the state sensitiviies @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fsrz(realtype t, int ie, N_Vector x, N_Vector *sx,
                     TempData *tdata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event timepoint, total derivative
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] sx Vector with the state sensitiviies @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fstau(realtype t, int ie, N_Vector x, N_Vector *sx,
                      TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Observables / measurements
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] user_data pointer to temp data object @type TempData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fy(realtype t, int it, N_Vector x, void *user_data,
                   ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of observables y w.r.t. model parameters p
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdydp(realtype t, int it, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of observables y w.r.t. state variables x
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdydx(realtype t, int it, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Event-resolved measurements
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fz(realtype t, int ie, N_Vector x, TempData *tdata,
                   ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of z, total derivative
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] sx Vector with the state sensitiviies @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fsz(realtype t, int ie, N_Vector x, N_Vector *sx,
                    TempData *tdata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurements z w.r.t. to model parameters
     * p
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdzdp(realtype t, int ie, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurements z w.r.t. to model states x
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdzdx(realtype t, int ie, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurements rz w.r.t. to model parameters
     * p
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdrzdp(realtype t, int ie, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurements rz w.r.t. to model states x
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in,out] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdrzdx(realtype t, int ie, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Right hand side of differential equation for states x
      * @param[in] t timepoint @type realtype
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] dx Vector with the derivative states (only DAE) @type
     * N_Vector
      * @param[out] xdot Vector with the right hand side @type N_Vector
      * @param[in] user_data pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                      void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
      * @param[in] user_data pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                       void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of dx/dt w.r.t. model parameters p
      * @param[in] t timepoint @type realtype
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] dx Vector with the derivative states (only DAE) @type
      * N_Vector
      * @param[in] user_data pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** State update functions for events
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in,out] x Vector with the states @type N_Vector
      * @param[in] xdot Vector with the right hand side @type N_Vector
      * @param[in] xdot_old Vector with the right hand side before the event
     * @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot,
                        N_Vector xdot_old, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity update functions for events, total derivative
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] xdot Vector with the right hand side @type N_Vector
      * @param[in] xdot_old Vector with the right hand side before the event
     * @type N_Vector
      * @param[in] sx Vector with the state sensitiviies @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot,
                         N_Vector xdot_old, N_Vector *sx, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Adjoint state update functions for events
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] xB Vector with the adjoint states @type N_Vector
      * @param[in] xdot Vector with the right hand side @type N_Vector
      * @param[in] xdot_old Vector with the right hand side before the event
     * @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB,
                         N_Vector xdot, N_Vector xdot_old, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Quadrature state update functions for events
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] xB Vector with the adjoint states @type N_Vector
      * @param[in] qBdot Vector with the adjoint quadrature states @type
     * N_Vector
      * @param[in] xdot Vector with the right hand side @type N_Vector
      * @param[in] xdot_old Vector with the right hand side before the event
     * @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB,
                         N_Vector qBdot, N_Vector xdot, N_Vector xdot_old,
                         TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Standard deviation of measurements
      * @param[in] t timepoint @type realtype
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fsigma_y(realtype t, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of standard deviation of measurements w.r.t. model
     * parameters p
      * @param[in] t timepoint @type realtype
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdsigma_ydp(realtype t, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Standard deviation of events
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fsigma_z(realtype t, int ie, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of standard deviation of events w.r.t. model parameters p
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] tdata pointer to temp data object @type TempData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdsigma_zdp(realtype t, int ie, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** negative log-likelihood of time-resolved measurements y
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fJy(realtype t, int it, N_Vector x, TempData *tdata,
                    const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** negative log-likelihood of event-resolved measurements z
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fJz(realtype t, int ie, N_Vector x, TempData *tdata,
                    const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** regularization of negative log-likelihood with roots of event-resolved
     * measurements rz
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fJrz(realtype t, int ie, N_Vector x, TempData *tdata,
                     const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. observables y
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJydy(realtype t, int it, N_Vector x, TempData *tdata,
                       const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. standard deviation sigma
      * @param[in] t timepoint @type realtype
      * @param[in] it timepoint index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJydsigma(realtype t, int it, N_Vector x, TempData *tdata,
                           const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurement negative log-likelihood Jz
     * w.r.t. event observables z
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJzdz(realtype t, int ie, N_Vector x, TempData *tdata,
                       const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigma
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJzdsigma(realtype t, int ie, N_Vector x, TempData *tdata,
                           const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurement negative log-likelihood
     * regularization Jrz w.r.t. event observables z
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJrzdz(realtype t, int ie, N_Vector x, TempData *tdata,
                        const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Sensitivity of event-resolved measurement negative log-likelihood
     * regularization Jrz w.r.t. standard deviation sigma
      * @param[in] t timepoint @type realtype
      * @param[in] ie event index @type int
      * @param[in] x Vector with the states @type N_Vector
      * @param[in] tdata pointer to temp data object @type TempData
      * @param[in] edata pointer to experimental data object @type ExpData
      * @param[in,out] rdata pointer to return data object @type ReturnData
      * @return status flag indicating successful execution @type int
      */
    virtual int fdJrzdsigma(realtype t, int ie, N_Vector x, TempData *tdata,
                            const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Second order Sensitivity of observables y w.r.t. state variables x
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in,out] tdata pointer to temp data object @type TempData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddydxdx(realtype t, int it, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /** Second order Sensitivity of observables y w.r.t. s.vs x and pars p
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in,out] tdata pointer to temp data object @type TempData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddydpdx(realtype t, int it, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Second order Sensitivity of observables y w.r.t. parameters p
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in,out] tdata pointer to temp data object @type TempData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddydpdp(realtype t, int it, N_Vector x, TempData *tdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Second order sensitivity of time-resolved measurement negative 
     * log-likelihood Jy w.r.t. observables y
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] tdata pointer to temp data object @type TempData
     * @param[in] edata pointer to experimental data object @type ExpData
     * @param[in,out] rdata pointer to return data object @type ReturnData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddJydydy(realtype t, int it, N_Vector x, TempData *tdata,
                       const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Second order sensitivity of time-resolved measurement negative
     * log-likelihood Jy w.r.t. standard deviation sigma and observables y
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] tdata pointer to temp data object @type TempData
     * @param[in] edata pointer to experimental data object @type ExpData
     * @param[in,out] rdata pointer to return data object @type ReturnData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddJydsigmady(realtype t, int it, N_Vector x, TempData *tdata,
                          const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Second order sensitivity of time-resolved measurement negative
     * log-likelihood Jy w.r.t. standard deviation sigma
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] tdata pointer to temp data object @type TempData
     * @param[in] edata pointer to experimental data object @type ExpData
     * @param[in,out] rdata pointer to return data object @type ReturnData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddJydsigmadsigma(realtype t, int it, N_Vector x, TempData *tdata,
                                  const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    /** Sensitivity of time-resolved measurement negative log-likelihood Jy 
     * w.r.t. standard deviation sigma in second order
     * @param[in] t timepoint @type realtype
     * @param[in] it timepoint index @type int
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] tdata pointer to temp data object @type TempData
     * @param[in] edata pointer to experimental data object @type ExpData
     * @param[in,out] rdata pointer to return data object @type ReturnData
     * @return status flag indicating successful execution @type int
     */
    virtual int fddJy_s2sigma(realtype t, int it, N_Vector x, TempData *tdata,
                                  const ExpData *edata, ReturnData *rdata) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    virtual int fdJdx(realtype t, N_Vector x, N_Vector dx, void *user_data){
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    virtual int fdJdp(realtype t, N_Vector x, N_Vector dx, void *user_data){
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }
    
    virtual int fddxdotdpdp(realtype t, N_Vector x, N_Vector dx, void *user_data){
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    /**
     * @brief Recurring terms in xdot
     * @param[in] t timepoint @type realtype
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states (only DAE) @type
     *N_Vector
     * @param[in] user_data object with model specifications @type TempData
     * @return status flag indicating successful execution @type int
     */
    virtual int fw(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
    virtual int fdwdp(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
    virtual int fdwdx(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
    virtual int fM(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
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
    virtual int fdfdx(realtype t, N_Vector x, N_Vector dx, void *user_data) {
        return AMICI_ERROR_NOT_IMPLEMENTED;
    }

    virtual ~Model();

    // Generic implementations
    int fsy(const int it, const TempData *tdata, ReturnData *rdata);

    int fsz_tf(const int ie, const TempData *tdata, ReturnData *rdata);

    int fsJy(const int it, const TempData *tdata, ReturnData *rdata);

    int fddJydpdp(const int it, TempData *tdata, const ExpData *edata,
               const ReturnData *rdata);
    
    int fdJydp(const int it, TempData *tdata, const ExpData *edata,
               const ReturnData *rdata);
    
    int fqBo2dot(realtype t, N_Vector x, N_Vector *sx, N_Vector xB,
                 N_Vector qBdot, void *user_data);

    int fdJydx(const int it, TempData *tdata, const ExpData *edata);

    int fsJz(const int ie, TempData *tdata, const ReturnData *rdata);

    int fdJzdp(const int ie, TempData *tdata, const ExpData *edata,
               const ReturnData *rdata);

    int fdJzdx(const int ie, TempData *tdata, const ExpData *edata);

    int initialize(const UserData *udata, TempData *tdata);

    int initializeStates(const double *x0data, TempData *tdata);

    int initHeaviside(TempData *tdata);

    /* Model dimensions */
    /** total number of model parameters */
    const int np;
    /** number of fixed parameters */
    const int nk;
    /** number of states */
    const int nx;
    /** number of states in the unaugmented system */
    const int nxtrue;
    /** number of observables */
    const int ny;
    /** number of observables in the unaugmented system */
    const int nytrue;
    /** number of event outputs */
    const int nz;
    /** number of event outputs in the unaugmented system */
    const int nztrue;
    /** number of events */
    const int ne;
    /** number of common expressions */
    const int nw;
    /** number of derivatives of common expressions wrt x */
    const int ndwdx;
    /** number of derivatives of common expressions wrt p */
    const int ndwdp;
    /** number of nonzero entries in jacobian */
    const int nnz;
    /** dimension of the augmented objective function for 2nd order ASA */
    const int nJ;
    /** upper bandwith of the jacobian */
    const int ubw;
    /** lower bandwith of the jacobian */
    const int lbw;
    /** flag indicating whether for sensi == AMICI_SENSI_ORDER_SECOND
     * directional or full second order derivative will be computed */
    const AMICI_o2mode o2mode;
    /** index indicating to which event an event output belongs */
    int *z2event = nullptr;
    /** flag array for DAE equations */
    realtype *idlist = nullptr;
};

#endif // MODEL_H
