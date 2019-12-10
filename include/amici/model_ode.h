#ifndef AMICI_MODEL_ODE_H
#define AMICI_MODEL_ODE_H

#include "amici/model.h"

#include <nvector/nvector_serial.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <utility>
#include <vector>

namespace amici {

class CVodeSolver;

/**
 * @brief The Model class represents an AMICI ODE model.
 *
 * The model does not contain any data, but represents the state
 * of the model at a specific time t. The states must not always be
 * in sync, but may be updated asynchroneously.
 */
class Model_ODE : public Model {
  public:
    /** default constructor */
    Model_ODE() = default;

    /**
     * @brief Constructor with model dimensions
     * @param nx_rdata number of state variables
     * @param nxtrue_rdata number of state variables of the non-augmented model
     * @param nx_solver number of state variables with conservation laws applied
     * @param nxtrue_solver number of state variables of the non-augmented model
     with conservation laws applied
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
     * @param ndxdotdw number of nonzero elements dxdotdw
     * @param ndJydy number of nonzero elements dJydy
     * @param nnz number of nonzero elements in Jacobian
     * @param ubw upper matrix bandwidth in the Jacobian
     * @param lbw lower matrix bandwidth in the Jacobian
     * @param o2mode second order sensitivity mode
     * @param p parameters
     * @param k constants
     * @param plist indexes wrt to which sensitivities are to be computed
     * @param idlist indexes indicating algebraic components (DAE only)
     * @param z2event mapping of event outputs to events
     * @param pythonGenerated flag indicating matlab or python wrapping
     * @param ndxdotdp_explicit number of nonzero elements dxdotdp_explicit
     * @param ndxdotdp_implicit number of nonzero elements dxdotdp_implicit
     */
    Model_ODE(const int nx_rdata, const int nxtrue_rdata, const int nx_solver,
              const int nxtrue_solver, const int ny, const int nytrue,
              const int nz, const int nztrue, const int ne, const int nJ,
              const int nw, const int ndwdx, const int ndwdp,
              const int ndxdotdw, std::vector<int> ndJydy,
              const int nnz, const int ubw, const int lbw,
              const SecondOrderMode o2mode, std::vector<realtype> const &p,
              std::vector<realtype> const &k, std::vector<int> const &plist,
              std::vector<realtype> const &idlist,
              std::vector<int> const &z2event, const bool pythonGenerated=false,
              const int ndxdotdp_explicit=0, const int ndxdotdp_implicit=0)
        : Model(nx_rdata, nxtrue_rdata, nx_solver, nxtrue_solver, ny, nytrue,
                nz, nztrue, ne, nJ, nw, ndwdx, ndwdp, ndxdotdw, std::move(ndJydy),
                nnz, ubw, lbw, o2mode, p, k, plist, idlist, z2event,
                pythonGenerated, ndxdotdp_explicit, ndxdotdp_implicit) {}

    void fJ(realtype t, realtype cj, const AmiVector &x, const AmiVector &dx,
            const AmiVector &xdot, SUNMatrix J) override;

    /**
     * @brief Implementation of fJ at the N_Vector level
     *
     * This function provides an
     * interface to the model specific routines for the solver
     * implementation as well as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     **/
    void fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J);

    /** implementation of fJB at the N_Vector level, this function provides an
     *interface to the model specific routines for the solver implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     **/
    void fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SUNMatrix JB);

    void fJSparse(realtype t, realtype cj, const AmiVector &x,
                  const AmiVector &dx, const AmiVector &xdot,
                  SUNMatrix J) override;

    /**
     * Implementation of fJSparse at the N_Vector level, this function
     * provides
     * an interface to the model specific routines for the solver implementation
     * aswell as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param J Matrix to which the Jacobian will be written
     */
    void fJSparse(realtype t, N_Vector x, SUNMatrix J);

    /** implementation of fJSparseB at the N_Vector level, this function
     * provides an interface to the model specific routines for the solver
     * implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     */
    void fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                   SUNMatrix JB);

    /** implementation of fJDiag at the N_Vector level, this function provides
     *an interface to the model specific routines for the solver implementation
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param x Vector with the states
     **/
    void fJDiag(realtype t, N_Vector JDiag, N_Vector x);

    /**
     * @brief diagonalized Jacobian (for preconditioning)
     * @param t timepoint
     * @param JDiag Vector to which the Jacobian diagonal will be written
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     **/
    void fJDiag(realtype t, AmiVector &JDiag, realtype cj, const AmiVector &x,
                const AmiVector &dx) override;

    void fJv(realtype t, const AmiVector &x, const AmiVector &dx,
             const AmiVector &xdot, const AmiVector &v, AmiVector &nJv,
             realtype cj) override;

    /** implementation of fJv at the N_Vector level.
     * @param t timepoint
     * @param x Vector with the states
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     * written
     **/
    void fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x);

    /**
     * @brief implementation of fJvB at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be written
     **/
    void fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB);

    void froot(realtype t, const AmiVector &x, const AmiVector &dx,
               gsl::span<realtype> root) override;

    /**
     * @brief implementation of froot at the N_Vector level
     *
     * This function provides an interface to the model specific routines for
     * the solver implementation aswell as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param root array with root function values
     */
    void froot(realtype t, N_Vector x, gsl::span<realtype> root);

    void fxdot(realtype t, const AmiVector &x, const AmiVector &dx,
               AmiVector &xdot) override;

    /** implementation of fxdot at the N_Vector level, this function provides an
     * interface to the model specific routines for the solver implementation
     * aswell as the AmiVector level implementation
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     */
    void fxdot(realtype t, N_Vector x, N_Vector xdot);

    /** implementation of fxBdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     */
    void fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot);

    /** implementation of fqBdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot);

    /** Sensitivity of dx/dt wrt model parameters w
     * @param t timepoint
     * @param x Vector with the states
     * @return status flag indicating successful execution
     */
    void fdxdotdw(realtype t, const N_Vector x);

    /** Explicit sensitivity of dx/dt wrt model parameters p
     * @param t timepoint
     * @param x Vector with the states
     * @return status flag indicating successful execution
     */
    void fdxdotdp(realtype t, const N_Vector x);

    void fdxdotdp(realtype t, const AmiVector &x, const AmiVector &dx) override;

    void fsxdot(realtype t, const AmiVector &x, const AmiVector &dx, int ip,
                const AmiVector &sx, const AmiVector &sdx,
                AmiVector &sxdot) override;

    /**
     * @brief implementation of fsxdot at the N_Vector level
     * @param t timepoint
     * @param x Vector with the states
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     */
    void fsxdot(realtype t, N_Vector x, int ip, N_Vector sx, N_Vector sxdot);

    std::unique_ptr<Solver> getSolver() override;

  protected:
    /** model specific implementation for fJ
     * @param J Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJ(realtype *J, realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w, const realtype *dwdx) = 0;

    /** model specific implementation for fJB
     * @param JB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJB(realtype *JB, realtype t, const realtype *x,
                     const realtype *p, const realtype *k, const realtype *h,
                     const realtype *xB, const realtype *w,
                     const realtype *dwdx);

    /** model specific implementation for fJSparse
     * @param JSparse Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparse(SUNMatrixContent_Sparse JSparse, realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h,
                          const realtype *w, const realtype *dwdx);

    /** model specific implementation for fJSparse, data only
     * @param JSparse Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparse(realtype *JSparse, realtype t, const realtype *x,
                          const realtype *p, const realtype *k,
                          const realtype *h, const realtype *w,
                          const realtype *dwdx);

    /**
     * @brief model specific implementation for fJSparse, column pointers
     * @param indexptrs column pointers
     **/
    virtual void fJSparse_colptrs(sunindextype *indexptrs);

    /**
     * @brief Model specific implementation for fJSparse, row values
     * @param indexvals row values
     **/
    virtual void fJSparse_rowvals(sunindextype *indexvals);

    /**
     * @brief Model specific implementation for fJSparseB
     * @param JSparseB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparseB(SUNMatrixContent_Sparse JSparseB, realtype t,
                           const realtype *x, const realtype *p,
                           const realtype *k, const realtype *h,
                           const realtype *xB, const realtype *w,
                           const realtype *dwdx);

    /** model specific implementation for fJSparseB
     * @param JSparseB data array
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparseB(realtype *JSparseB, realtype t, const realtype *x,
                           const realtype *p, const realtype *k,
                           const realtype *h, const realtype *xB,
                           const realtype *w, const realtype *dwdx);

    /**
     * @brief Model specific implementation for fJSparse, column pointers
     * @param indexptrs column pointers
     **/
    virtual void fJSparseB_colptrs(sunindextype *indexptrs);

    /**
     * @brief Model specific implementation for fJSparse, row values
     * @param indexvals row values
     **/
    virtual void fJSparseB_rowvals(sunindextype *indexvals);

    /**
     * @brief Model specific implementation for fJDiag
     * @param JDiag Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJDiag(realtype *JDiag, realtype t, const realtype *x,
                        const realtype *p, const realtype *k, const realtype *h,
                        const realtype *w, const realtype *dwdx);

    /**
     * @brief model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     **/
    virtual void froot(realtype *root, realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h);

    /** model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     **/
    virtual void fxdot(realtype *xdot, realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w) = 0;

    /** model specific implementation of fdxdotdp, with w chainrule (Matlab)
     * @param dxdotdp partial derivative xdot wrt p
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param ip parameter index
     * @param w vector with helper variables
     * @param dwdp derivative of w wrt p
     */
    virtual void fdxdotdp(realtype *dxdotdp, realtype t, const realtype *x,
                          const realtype *p, const realtype *k,
                          const realtype *h, int ip, const realtype *w,
                          const realtype *dwdp);

    /** model specific implementation of fdxdotdp_explicit, no w chainrule (Py)
     * @param dxdotdp_explicit partial derivative xdot wrt p
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     */
    virtual void fdxdotdp_explicit(realtype *dxdotdp_explicit, realtype t,
                                   const realtype *x, const realtype *p,
                                   const realtype *k, const realtype *h,
                                   const realtype *w);

    /** model specific implementation of fdxdotdp_explicit, colptrs part
     * @param indexptrs column pointers
     */
    virtual void fdxdotdp_explicit_colptrs(sunindextype *indexptrs);

    /** model specific implementation of fdxdotdp_explicit, rowvals part
     * @param indexvals row values
     */
    virtual void fdxdotdp_explicit_rowvals(sunindextype *indexvals);

    /** model specific implementation of fdxdotdp_implicit, colptrs part
     * @param indexptrs column pointers
     */
    virtual void fdxdotdp_implicit_colptrs(sunindextype *indexptrs);

    /** model specific implementation of fdxdotdp_implicit, rowvals part
     * @param indexvals row values
     */
    virtual void fdxdotdp_implicit_rowvals(sunindextype *indexvals);

    /** model specific implementation of fdxdotdw, data part
     * @param dxdotdw partial derivative xdot wrt w
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     */
    virtual void fdxdotdw(realtype *dxdotdw, realtype t, const realtype *x,
                          const realtype *p, const realtype *k,
                          const realtype *h, const realtype *w);

    /** model specific implementation of fdxdotdw, colptrs part
     * @param indexptrs column pointers
     */
    virtual void fdxdotdw_colptrs(sunindextype *indexptrs);

    /** model specific implementation of fdxdotdw, rowvals part
     * @param indexvals row values
     */
    virtual void fdxdotdw_rowvals(sunindextype *indexvals);
};
} // namespace amici

#endif // MODEL_H
