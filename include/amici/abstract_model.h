#ifndef AMICI_ABSTRACT_MODEL_H
#define AMICI_ABSTRACT_MODEL_H

#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <memory>

namespace amici {

class Solver;

/**
 * @brief Abstract base class of amici::Model defining functions that need to
 * be implemented in an AMICI model.
 *
 * Some functions have empty default implementations or throw.
 * This class shall not have any data members.
 */
class AbstractModel {
  public:

    virtual ~AbstractModel() = default;

    /**
     * Retrieves the solver object
     * @return The Solver instance
     */
    virtual std::unique_ptr<Solver> getSolver() = 0;

    /**
     * Root function
     * @param t time
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param root array to which values of the root function will be written
     */
    virtual void froot(realtype t, AmiVector *x, AmiVector *dx,
                       realtype *root) = 0;

    /**
     * Residual function
     * @param t time
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param xdot array to which values of the residual function will be
     * written
     */
    virtual void fxdot(realtype t, AmiVector *x, AmiVector *dx,
                       AmiVector *xdot) = 0;

    /**
     * Sensitivity Residual function
     * @param t time
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param ip parameter index
     * @param sx sensitivity state
     * @param sdx time derivative of sensitivity state (DAE only)
     * @param sxdot array to which values of the sensitivity residual function
     * will be written
     */
    virtual void fsxdot(realtype t, AmiVector *x, AmiVector *dx, int ip,
                        AmiVector *sx, AmiVector *sdx, AmiVector *sxdot) = 0;

    /**
     * Dense Jacobian function
     * @param t time
     * @param cj scaling factor (inverse of timestep, DAE only)
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param xdot values of residual function (unused)
     * @param J dense matrix to which values of the jacobian will be written
     */
    virtual void fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                    AmiVector *xdot, DlsMat J) = 0;

    /**
     * Sparse Jacobian function
     * @param t time
     * @param cj scaling factor (inverse of timestep, DAE only)
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param xdot values of residual function (unused)
     * @param J sparse matrix to which values of the Jacobian will be written
     */
    virtual void fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                          AmiVector *xdot, SlsMat J) = 0;

    /**
     * Diagonal Jacobian function
     * @param t time
     * @param Jdiag array to which the diagonal of the Jacobian will be written
     * @param cj scaling factor (inverse of timestep, DAE only)
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @return flag indicating successful evaluation
     */
    virtual void fJDiag(realtype t, AmiVector *Jdiag, realtype cj, AmiVector *x,
                        AmiVector *dx) = 0;

    /**
     * Parameter derivative of residual function
     * @param t time
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @return flag indicating successful evaluation
     */
    virtual void fdxdotdp(realtype t, AmiVector *x, AmiVector *dx) = 0;

    /**
     * Jacobian multiply function
     * @param t time
     * @param x state
     * @param dx time derivative of state (DAE only)
     * @param xdot values of residual function (unused)
     * @param v multiplication vector (unused)
     * @param nJv array to which result of multiplication will be written
     * @param cj scaling factor (inverse of timestep, DAE only)
     */
    virtual void fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                     AmiVector *v, AmiVector *nJv, realtype cj) = 0;

    /**
     * @brief Returns the amici version that was used to generate the model
     * @return ver amici version string
     */
    virtual const std::string getAmiciVersion() const;

    /**
     * @brief Returns the amici commit that was used to generate the model
     * @return ver amici commit string
     */
    virtual const std::string getAmiciCommit() const;

    /**
     * Model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p,
                     const realtype *k);

    /**
     * Function indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     * @return flag inidication whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    virtual bool isFixedParameterStateReinitializationAllowed() const;

    /**
     * Model specific implementation of fx0_fixedParameters
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0_fixedParameters(realtype *x0, const realtype t,
                                     const realtype *p, const realtype *k);

    /**
     * Model specific implementation of fsx0_fixedParameters
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0_fixedParameters(realtype *sx0, const realtype t,
                                      const realtype *x0, const realtype *p,
                                      const realtype *k, const int ip);

    /**
     * Model specific implementation of fsx0
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0(realtype *sx0, const realtype t, const realtype *x0,
                      const realtype *p, const realtype *k, const int ip);

    /**
     * Initial value for time derivative of states (only necessary for DAEs)
     * @param x0 Vector with the initial states
     * @param dx0 Vector to which the initial derivative states will be
     * written (only DAE)
     **/
    virtual void fdx0(AmiVector *x0, AmiVector *dx0);

    /**
     * Model specific implementation of fstau
     * @param stau total derivative of event timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     * @param ie event index
     **/
    virtual void fstau(realtype *stau, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *sx, const int ip, const int ie);

    /**
     * Model specific implementation of fy
     * @param y model output at current timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     **/
    virtual void fy(realtype *y, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *w);

    /**
     * Model specific implementation of fdydp
     * @param dydp partial derivative of observables y w.r.t. model parameters p
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     * @param w repeating elements vector
     * @param dwdp Recurring terms in xdot, parameter derivative
     **/
    virtual void fdydp(realtype *dydp, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const int ip, const realtype *w, const realtype *dwdp);

    /**
     * Model specific implementation of fdydx
     * @param dydx partial derivative of observables y w.r.t. model states x
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     * @param dwdx Recurring terms in xdot, state derivative
     **/
    virtual void fdydx(realtype *dydx, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *dwdx);

    /**
     * Model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fz(realtype *z, const int ie, const realtype t,
                    const realtype *x, const realtype *p, const realtype *k,
                    const realtype *h);

    /**
     * Model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     **/
    virtual void fsz(realtype *sz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h, const realtype *sx, const int ip);

    /**
     * Model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     *not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void frz(realtype *rz, const int ie, const realtype t,
                     const realtype *x, const realtype *p, const realtype *k,
                     const realtype *h);

    /**
     * Model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heavyside vector
     * @param ip sensitivity index
     **/
    virtual void fsrz(realtype *srz, const int ie, const realtype t,
                      const realtype *x, const realtype *p, const realtype *k,
                      const realtype *h, const realtype *sx, const int ip);

    /**
     * Model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     *parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdzdp(realtype *dzdp, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h, const int ip);

    /**
     * Model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     *states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdzdx(realtype *dzdx, const int ie, const realtype t,
                       const realtype *x, const realtype *p, const realtype *k,
                       const realtype *h);

    /**
     * Model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     *p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h, const int ip);

    /**
     * Model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                        const realtype *x, const realtype *p, const realtype *k,
                        const realtype *h);

    /**
     * Model specific implementation of fdeltax
     * @param deltax state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     **/
    virtual void fdeltax(realtype *deltax, const realtype t, const realtype *x,
                         const realtype *p, const realtype *k,
                         const realtype *h, const int ie, const realtype *xdot,
                         const realtype *xdot_old);

    /**
     * Model specific implementation of fdeltasx
     * @param deltasx sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param sx state sensitivity
     * @param stau event-time sensitivity
     **/
    virtual void fdeltasx(realtype *deltasx, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h,
                          const realtype *w, const int ip, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *sx, const realtype *stau);

    /**
     * Model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     **/
    virtual void fdeltaxB(realtype *deltaxB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ie,
                          const realtype *xdot, const realtype *xdot_old,
                          const realtype *xB);

    /**
     * Model specific implementation of fdeltaqB
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     **/
    virtual void fdeltaqB(realtype *deltaqB, const realtype t,
                          const realtype *x, const realtype *p,
                          const realtype *k, const realtype *h, const int ip,
                          const int ie, const realtype *xdot,
                          const realtype *xdot_old, const realtype *xB);

    /**
     * Model specific implementation of fsigmay
     * @param sigmay standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p,
                         const realtype *k);

    /**
     * Model specific implementation of fsigmay
     * @param dsigmaydp partial derivative of standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t,
                            const realtype *p, const realtype *k, const int ip);

    /**
     * Model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                         const realtype *k);

    /**
     * Model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     *measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                            const realtype *p, const realtype *k, const int ip);

    /**
     * Model specific implementation of fJy
     * @param nllh negative log-likelihood for measurements y
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurements at timepoint
     **/
    virtual void fJy(realtype *nllh, const int iy, const realtype *p,
                     const realtype *k, const realtype *y,
                     const realtype *sigmay, const realtype *my);

    /**
     * Model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     **/
    virtual void fJz(realtype *nllh, const int iz, const realtype *p,
                     const realtype *k, const realtype *z,
                     const realtype *sigmaz, const realtype *mz);

    /**
     * Model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fJrz(realtype *nllh, const int iz, const realtype *p,
                      const realtype *k, const realtype *z,
                      const realtype *sigmaz);

    /**
     * Model specific implementation of fdJydy
     * @param dJydy partial derivative of time-resolved measurement negative
     *log-likelihood Jy
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydy(realtype *dJydy, const int iy, const realtype *p,
                        const realtype *k, const realtype *y,
                        const realtype *sigmay, const realtype *my);

    /**
     * Model specific implementation of fdJydsigma
     * @param dJydsigma Sensitivity of time-resolved measurement
     * negative log-likelihood Jy w.r.t. standard deviation sigmay
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydsigma(realtype *dJydsigma, const int iy,
                            const realtype *p, const realtype *k,
                            const realtype *y, const realtype *sigmay,
                            const realtype *my);

    /**
     * Model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     *log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                        const realtype *k, const realtype *z,
                        const realtype *sigmaz, const realtype *mz);

    /**
     * Model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdsigma(realtype *dJzdsigma, const int iz,
                            const realtype *p, const realtype *k,
                            const realtype *z, const realtype *sigmaz,
                            const realtype *mz);

    /**
     * Model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                         const realtype *k, const realtype *rz,
                         const realtype *sigmaz);

    /**
     * Model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                             const realtype *p, const realtype *k,
                             const realtype *rz, const realtype *sigmaz);

    /**
     * Model specific implementation of fw
     * @param w Recurring terms in xdot
     * @param t timepoint
     * @param x vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param tcl total abundances for conservations laws
     */
    virtual void fw(realtype *w, const realtype t, const realtype *x,
                    const realtype *p, const realtype *k, const realtype *h,
                    const realtype *tcl);

    /**
     * Model specific implementation of dwdp
     * @param dwdp Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param tcl total abundances for conservations laws
     * @param stcl sensitivities of total abundances for conservations laws
     */
    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *tcl,
                       const realtype *stcl);

    /**
     * Model specific implementation of dwdx
     * @param dwdx Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param tcl total abundances for conservations laws
     */
    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h,
                       const realtype *w, const realtype *tcl);
};

} // namespace amici

#endif // AMICI_ABSTRACT_MODEL_H
