#ifndef AMICI_MODEL_DAE_H
#define AMICI_MODEL_DAE_H

#include "amici/model.h"

#include <nvector/nvector_serial.h>

#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <utility>
#include <vector>

namespace amici {

class ExpData;
class IDASolver;

/**
 * @brief The Model class represents an AMICI DAE model.
 *
 * The model does not contain any data, but represents the state
 * of the model at a specific time t. The states must not always be
 * in sync, but may be updated asynchronously.
 */
class Model_DAE : public Model {
  public:
    /** default constructor */
    Model_DAE() = default;

    /**
     * @brief Constructor with model dimensions
     * @param model_dimensions Model dimensions
     * @param simulation_parameters Simulation parameters
     * @param o2mode second order sensitivity mode
     * @param idlist indexes indicating algebraic components (DAE only)
     * @param z2event mapping of event outputs to events
     * @param pythonGenerated flag indicating matlab or python wrapping
     * @param ndxdotdp_explicit number of nonzero elements dxdotdp_explicit
     */
    Model_DAE(const ModelDimensions &model_dimensions,
              SimulationParameters simulation_parameters,
              const SecondOrderMode o2mode,
              std::vector<realtype> const &idlist,
              std::vector<int> const &z2event, const bool pythonGenerated=false,
              const int ndxdotdp_explicit=0)
        : Model(model_dimensions, simulation_parameters,
                o2mode, idlist, z2event, pythonGenerated,
                ndxdotdp_explicit) {
        derived_state_.M_ = SUNMatrixWrapper(nx_solver, nx_solver);
    }

    void fJ(realtype t, realtype cj, const AmiVector &x, const AmiVector &dx,
            const AmiVector &xdot, SUNMatrix J) override;

    /**
     * @brief Jacobian of xdot with respect to states x
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     **/
    void fJ(realtype t, realtype cj, const_N_Vector x, const_N_Vector dx,
            const_N_Vector xdot, SUNMatrix J);

    void fJB(const realtype t, realtype cj, const AmiVector &x,
             const AmiVector &dx, const AmiVector &xB, const AmiVector &dxB,
             const AmiVector &xBdot, SUNMatrix JB) override;

    /**
     * @brief Jacobian of xBdot with respect to adjoint state xB
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param JB Matrix to which the Jacobian will be written
     **/
    void fJB(realtype t, realtype cj, const_N_Vector x, const_N_Vector dx,
             const_N_Vector xB, const_N_Vector dxB, SUNMatrix JB);

    void fJSparse(realtype t, realtype cj, const AmiVector &x,
                  const AmiVector &dx, const AmiVector &xdot,
                  SUNMatrix J) override;

    /**
     * @brief J in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param J Matrix to which the Jacobian will be written
     */
    void fJSparse(realtype t, realtype cj, const_N_Vector x, const_N_Vector dx,
                  SUNMatrix J);

    void fJSparseB(const realtype t, realtype cj, const AmiVector &x,
                   const AmiVector &dx, const AmiVector &xB,
                   const AmiVector &dxB, const AmiVector &xBdot,
                   SUNMatrix JB) override;

    /**
     * @brief JB in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param JB Matrix to which the Jacobian will be written
     */
    void fJSparseB(realtype t, realtype cj, const_N_Vector x, const_N_Vector dx,
                   const_N_Vector xB, const_N_Vector dxB, SUNMatrix JB);

    /**
     * @brief Diagonal of the Jacobian (for preconditioning)
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

    /**
     * @brief Matrix vector product of J with a vector v (for iterative solvers)
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     * written
     **/
    void fJv(realtype t, const_N_Vector x, const_N_Vector dx, const_N_Vector v,
             N_Vector Jv, realtype cj);

    /**
     * @brief Matrix vector product of JB with a vector v (for iterative solvers)
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be written
     * @param cj scalar in Jacobian (inverse stepsize)
     **/

    void fJvB(realtype t, const_N_Vector x, const_N_Vector dx,
              const_N_Vector xB, const_N_Vector dxB,
              const_N_Vector vB, N_Vector JvB, realtype cj);

    void froot(realtype t, const AmiVector &x, const AmiVector &dx,
               gsl::span<realtype> root) override;

    /**
     * @brief Event trigger function for events
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param root array with root function values
     */
    void froot(realtype t, const_N_Vector x, const_N_Vector dx, gsl::span<realtype> root);

    void fxdot(realtype t, const AmiVector &x, const AmiVector &dx,
               AmiVector &xdot) override;

    /**
     * @brief Residual function of the DAE
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     */
    void fxdot(realtype t, const_N_Vector x, const_N_Vector dx, N_Vector xdot);

    /**
     * @brief Right hand side of differential equation for adjoint state xB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     */
    void fxBdot(realtype t, const_N_Vector x, const_N_Vector dx,
                const_N_Vector xB, const_N_Vector dxB, N_Vector xBdot);

    /**
     * @brief Right hand side of integral equation for quadrature states qB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void fqBdot(realtype t, const_N_Vector x, const_N_Vector dx,
                const_N_Vector xB, const_N_Vector dxB,
                N_Vector qBdot);

    void fxBdot_ss(const realtype t, const AmiVector &xB,
                   const AmiVector &dxB, AmiVector &xBdot) override;

    /**
     * @brief Implementation of fxBdot for steady state case at the N_Vector level
     * @param t timepoint
     * @param xB Vector with the adjoint state
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     */
    void fxBdot_ss(realtype t, const_N_Vector xB, const_N_Vector dxB,
                   N_Vector xBdot) const;

    /**
     * @brief Implementation of fqBdot for steady state at the N_Vector level
     * @param t timepoint
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param qBdot Vector with the adjoint quadrature right hand side
     */
    void fqBdot_ss(realtype t, const_N_Vector xB, const_N_Vector dxB,
                   N_Vector qBdot) const;

    /**
     * @brief Sparse Jacobian function backward, steady state case
     * @param JB sparse matrix to which values of the Jacobian will be written
     */
    void fJSparseB_ss(SUNMatrix JB) override;

    /**
     * @brief Computes the sparse backward Jacobian for steadystate integration
     * and writes it to the model member
     * @param t timepoint
     * @param cj scalar in Jacobian
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint state right hand side
     */
    void writeSteadystateJB(const realtype t, realtype cj,
                            const AmiVector &x, const AmiVector &dx,
                            const AmiVector &xB, const AmiVector &dxB,
                            const AmiVector &xBdot) override;

    /**
     * @brief Sensitivity of dx/dt wrt model parameters p
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     */
    void fdxdotdp(realtype t, const const_N_Vector x, const const_N_Vector dx);
    void fdxdotdp(const realtype t, const AmiVector &x,
                  const AmiVector &dx) override {
        fdxdotdp(t, x.getNVector(), dx.getNVector());
    };

    void fsxdot(realtype t, const AmiVector &x, const AmiVector &dx, int ip,
                const AmiVector &sx, const AmiVector &sdx,
                AmiVector &sxdot) override;
    /**
     * @brief Right hand side of differential equation for state sensitivities sx
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sdx Vector with the derivative state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     */
    void fsxdot(realtype t, const_N_Vector x, const_N_Vector dx, int ip,
                const_N_Vector sx, const_N_Vector sdx, N_Vector sxdot);

    /**
     * @brief Mass matrix for DAE systems
     * @param t timepoint
     * @param x Vector with the states
     */
    void fM(realtype t, const_N_Vector x);

    std::unique_ptr<Solver> getSolver() override;

  protected:
    /**
     * @brief Model specific implementation for fJSparse
     * @param JSparse Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h Heaviside vector
     * @param cj scaling factor, inverse of the step size
     * @param dx Vector with the derivative states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparse(SUNMatrixContent_Sparse JSparse, realtype t,
                          const realtype *x, const double *p, const double *k,
                          const realtype *h, realtype cj, const realtype *dx,
                          const realtype *w, const realtype *dwdx) = 0;

    /**
     * @brief Model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h Heaviside vector
     * @param dx Vector with the derivative states
     **/
    virtual void froot(realtype *root, realtype t, const realtype *x,
                       const double *p, const double *k, const realtype *h,
                       const realtype *dx);

    /**
     * @brief Model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h Heaviside vector
     * @param w vector with helper variables
     * @param dx Vector with the derivative states
     **/
    virtual void fxdot(realtype *xdot, realtype t, const realtype *x,
                       const double *p, const double *k, const realtype *h,
                       const realtype *dx, const realtype *w) = 0;

    /**
     * @brief Model specific implementation of fdxdotdp
     * @param dxdotdp partial derivative xdot wrt p
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h Heaviside vector
     * @param ip parameter index
     * @param dx Vector with the derivative states
     * @param w vector with helper variables
     * @param dwdp derivative of w wrt p
     */
    virtual void fdxdotdp(realtype *dxdotdp, realtype t, const realtype *x,
                          const realtype *p, const realtype *k,
                          const realtype *h, int ip, const realtype *dx,
                          const realtype *w, const realtype *dwdp);

    /**
     * @brief Model specific implementation of fM
     * @param M mass matrix
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     */
    virtual void fM(realtype *M, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k);
};
} // namespace amici

#endif // MODEL_H
