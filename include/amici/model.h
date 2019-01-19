#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

#include "amici/abstract_model.h"
#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"

#include <memory>
#include <vector>

namespace amici {

class ReturnData;
class ExpData;
class Model;
class Solver;

} // namespace amici

// for serialization friend in amici::Model
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::Model &u, const unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief The Model class represents an AMICI ODE model.
 * The model can compute various model related quantities based
 * on symbolically generated code.
 */
class Model : public AbstractModel {
  public:
    /** default constructor */
    Model();

    /**
     * Constructor with model dimensions
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
     * @param nnz number of nonzero elements in Jacobian
     * @param ubw upper matrix bandwidth in the Jacobian
     * @param lbw lower matrix bandwidth in the Jacobian
     * @param o2mode second order sensitivity mode
     * @param p parameters
     * @param k constants
     * @param plist indexes wrt to which sensitivities are to be computed
     * @param idlist indexes indicating algebraic components (DAE only)
     * @param z2event mapping of event outputs to events
     */
    Model(const int nx_rdata, const int nxtrue_rdata, const int nx_solver,
          const int nxtrue_solver, const int ny, const int nytrue, const int nz,
          const int nztrue, const int ne, const int nJ, const int nw,
          const int ndwdx, const int ndwdp, const int nnz, const int ubw,
          const int lbw, amici::SecondOrderMode o2mode,
          const std::vector<amici::realtype> &p, std::vector<amici::realtype> k,
          const std::vector<int> &plist, std::vector<amici::realtype> idlist,
          std::vector<int> z2event);

    /** destructor */
    virtual ~Model() = default;

    /**
     * Copy assignment is disabled until const members are removed
     * @param other object to copy from
     * @return
     */
    Model &operator=(Model const &other) = delete;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Model *clone() const = 0;

    // Overloaded base class methods
    using AbstractModel::fdeltaqB;
    using AbstractModel::fdeltasx;
    using AbstractModel::fdeltax;
    using AbstractModel::fdeltaxB;
    using AbstractModel::fdJrzdsigma;
    using AbstractModel::fdJrzdz;
    using AbstractModel::fdJydsigma;
    using AbstractModel::fdJydy;
    using AbstractModel::fdJzdsigma;
    using AbstractModel::fdJzdz;
    using AbstractModel::fdrzdp;
    using AbstractModel::fdrzdx;
    using AbstractModel::fdsigmaydp;
    using AbstractModel::fdsigmazdp;
    using AbstractModel::fdwdp;
    using AbstractModel::fdwdx;
    using AbstractModel::fdydp;
    using AbstractModel::fdydx;
    using AbstractModel::fdzdp;
    using AbstractModel::fdzdx;
    using AbstractModel::fJrz;
    using AbstractModel::fJy;
    using AbstractModel::fJz;
    using AbstractModel::frz;
    using AbstractModel::fsigmay;
    using AbstractModel::fsigmaz;
    using AbstractModel::fsrz;
    using AbstractModel::fstau;
    using AbstractModel::fsx0;
    using AbstractModel::fsx0_fixedParameters;
    using AbstractModel::fsz;
    using AbstractModel::fw;
    using AbstractModel::fx0;
    using AbstractModel::fx0_fixedParameters;
    using AbstractModel::fy;
    using AbstractModel::fz;

    /**
     * Expands conservation law for states
     * @param x_rdata pointer to state variables with conservation laws
     * expanded (stored in rdata)
     * @param x_solver pointer to state variables with conservation laws
     * applied (solver returns this)
     */
    void fx_rdata(AmiVector *x_rdata, const AmiVector *x_solver);

    /**
     * Expands conservation law for state sensitivities
     * @param sx_rdata pointer to state variable sensitivities with
     * conservation laws expanded (stored in rdata)
     * @param sx_solver pointer to state variable sensitivities with
     * conservation laws applied (solver returns this)
     */
    void fsx_rdata(AmiVectorArray *sx_rdata, const AmiVectorArray *sx_solver);

    /**
     * Initial states
     * @param x pointer to state variables
     */
    void fx0(AmiVector *x);

    /**
     * Sets only those initial states that are specified via fixedParmeters
     * @param x pointer to state variables
     */
    void fx0_fixedParameters(AmiVector *x);

    /**
     * Initial value for initial state sensitivities
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     **/
    void fsx0(AmiVectorArray *sx, const AmiVector *x);

    /**
     * Sets only those initial states sensitivities that are affected from fx0
     *fixedParmeters
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     **/
    void fsx0_fixedParameters(AmiVectorArray *sx, const AmiVector *x);

    /**
     * Sensitivity of derivative initial states sensitivities sdx0 (only
     *  necessary for DAEs)
     **/
    virtual void fsdx0();

    /**
     * Sensitivity of event timepoint, total derivative
     * @param t current timepoint
     * @param ie event index
     * @param x pointer to state variables
     * @param sx pointer to state sensitivity variables
     */
    void fstau(const realtype t, const int ie, const AmiVector *x,
               const AmiVectorArray *sx);

    /**
     * Observables / measurements
     * @param t current timepoint
     * @param it timepoint index
     * @param x current state
     * @param rdata pointer to return data instance
     */
    void fy(const realtype t, const int it, const AmiVector *x,
            ReturnData *rdata);

    /**
     * Partial derivative of observables y w.r.t. model parameters p
     * @param t current timepoint
     * @param x current state
     */
    void fdydp(const realtype t, const AmiVector *x);

    /**
     * Partial derivative of observables y w.r.t. state variables x
     * @param t current timepoint
     * @param x current state
     */
    void fdydx(const realtype t, const AmiVector *x);

    /** Event-resolved output
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
    void fz(const int nroots, const int ie, const realtype t,
            const AmiVector *x, ReturnData *rdata);

    /** Sensitivity of z, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
    void fsz(const int nroots, const int ie, const realtype t,
             const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);

    /**
     * Event root function of events (equal to froot but does not include
     * non-output events)
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param rdata pointer to return data instance
     */
    void frz(const int nroots, const int ie, const realtype t,
             const AmiVector *x, ReturnData *rdata);

    /**
     * Sensitivity of rz, total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivities
     * @param rdata pointer to return data instance
     */
    void fsrz(const int nroots, const int ie, const realtype t,
              const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);

    /**
     * Partial derivative of event-resolved output z w.r.t. to model parameters
     * p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdzdp(const realtype t, const int ie, const AmiVector *x);

    /**
     * Partial derivative of event-resolved output z w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdzdx(const realtype t, const int ie, const AmiVector *x);

    /**
     * Sensitivity of event-resolved root output w.r.t. to model parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdrzdp(const realtype t, const int ie, const AmiVector *x);

    /**
     * Sensitivity of event-resolved measurements rz w.r.t. to model states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdrzdx(const realtype t, const int ie, const AmiVector *x);

    /**
     * State update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void fdeltax(const int ie, const realtype t, const AmiVector *x,
                 const AmiVector *xdot, const AmiVector *xdot_old);

    /**
     * Sensitivity update functions for events, total derivative
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param sx current state sensitivity
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void fdeltasx(const int ie, const realtype t, const AmiVector *x,
                  const AmiVectorArray *sx, const AmiVector *xdot,
                  const AmiVector *xdot_old);

    /**
     * Adjoint state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void fdeltaxB(const int ie, const realtype t, const AmiVector *x,
                  const AmiVector *xB, const AmiVector *xdot,
                  const AmiVector *xdot_old);

    /**
     * Quadrature state update functions for events
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void fdeltaqB(const int ie, const realtype t, const AmiVector *x,
                  const AmiVector *xB, const AmiVector *xdot,
                  const AmiVector *xdot_old);

    /**
     * Standard deviation of measurements
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void fsigmay(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Partial derivative of standard deviation of measurements w.r.t. model
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to ExpData data instance holding sigma values
     */
    void fdsigmaydp(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Standard deviation of events
     * @param t current timepoint
     * @param ie event index
     * @param nroots array with event numbers
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void fsigmaz(const realtype t, const int ie, const int *nroots,
                 ReturnData *rdata, const ExpData *edata);

    /**
     * Sensitivity of standard deviation of events measurements w.r.t. model
     * parameters p
     * @param t current timepoint
     * @param ie event index
     * @param nroots array with event numbers
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdsigmazdp(const realtype t, const int ie, const int *nroots,
                    ReturnData *rdata, const ExpData *edata);

    /**
     * Negative log-likelihood of measurements y
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fJy(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Negative log-likelihood of event-resolved measurements z
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fJz(const int nroots, ReturnData *rdata, const ExpData *edata);

    /**
     * Regularization of negative log-likelihood with roots of event-resolved
     * measurements rz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fJrz(const int nroots, ReturnData *rdata, const ExpData *edata);

    /**
     * Partial derivative of time-resolved measurement negative log-likelihood
     * Jy
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJydy(const int it, const ReturnData *rdata, const ExpData *edata);

    /**
     * Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. standard deviation sigma
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJydsigma(const int it, const ReturnData *rdata,
                    const ExpData *edata);

    /**
     * Partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJzdz(const int nroots, const ReturnData *rdata,
                const ExpData *edata);

    /**
     * Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJzdsigma(const int nroots, const ReturnData *rdata,
                    const ExpData *edata);

    /**
     * Partial derivative of event measurement negative log-likelihood Jz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJrzdz(const int nroots, const ReturnData *rdata,
                 const ExpData *edata);

    /**
     * Sensitivity of event measurement negative log-likelihood Jz
     * w.r.t. standard deviation sigmaz
     * @param nroots event index
     * @param rdata pointer to return data instance
     * @param edata pointer to experimental data instance
     */
    void fdJrzdsigma(const int nroots, const ReturnData *rdata,
                     const ExpData *edata);

    /**
     * Sensitivity of measurements y, total derivative sy = dydx * sx + dydp
     * @param it timepoint index
     * @param sx pointer to state sensitivities
     * @param rdata pointer to return data instance
     */
    void fsy(const int it, const AmiVectorArray *sx, ReturnData *rdata);

    /**
     * Sensitivity of z at final timepoint (ignores sensitivity of timepoint),
     * total derivative
     * @param nroots number of events for event index
     * @param ie event index
     * @param rdata pointer to return data instance
     */
    void fsz_tf(const int *nroots, const int ie, ReturnData *rdata);

    /**
     * Sensitivity of time-resolved measurement negative log-likelihood Jy,
     * total derivative
     * @param it timepoint index
     * @param sx pointer to state sensitivities
     * @param dJydx vector with values of state derivative of Jy
     * @param rdata pointer to return data instance
     */
    void fsJy(const int it, const std::vector<realtype> &dJydx,
              const AmiVectorArray *sx, ReturnData *rdata);

    /**
     * Compute sensitivity of time-resolved measurement negative log-likelihood
     * Jy w.r.t. parameters for the given timepoint. Add result to respective
     * fields in rdata.
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void fdJydp(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Sensitivity of time-resolved measurement negative log-likelihood Jy
     * w.r.t. state variables
     * @param dJydx pointer to vector with values of state derivative of Jy
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
    void fdJydx(std::vector<realtype> *dJydx, const int it,
                const ExpData *edata);

    /**
     * Sensitivity of event-resolved measurement negative log-likelihood Jz,
     * total derivative
     * @param nroots event index
     * @param dJzdx vector with values of state derivative of Jz
     * @param sx pointer to state sensitivities
     * @param rdata pointer to return data instance
     */
    void fsJz(const int nroots, const std::vector<realtype> &dJzdx,
              const AmiVectorArray *sx, ReturnData *rdata);

    /**
     * Sensitivity of event-resolved measurement negative log-likelihood Jz
     * w.r.t. parameters
     * @param nroots event index
     * @param t current timepoint
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void fdJzdp(const int nroots, realtype t, const ExpData *edata,
                const ReturnData *rdata);

    /**
     * Sensitivity of event-resolved measurement negative log-likelihood Jz
     * w.r.t. state variables
     * @param dJzdx pointer to vector with values of state derivative of Jz
     * @param nroots event index
     * @param t current timepoint
     * @param edata pointer to experimental data instance
     * @param rdata pointer to return data instance
     */
    void fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t,
                const ExpData *edata, const ReturnData *rdata);

    /**
     * Initialization of model properties
     * @param x pointer to state variables
     * @param dx pointer to time derivative of states (DAE only)
     * @param sx pointer to state variable sensititivies
     * @param sdx pointer to time derivative of state sensitivities
     * (DAE only)
     * @param computeSensitivities flag indicating whether sensitivities
     * are to be computed
     */
    void initialize(AmiVector *x, AmiVector *dx, AmiVectorArray *sx,
                    AmiVectorArray *sdx, bool computeSensitivities);

    /**
     * Initialization of initial states
     * @param x pointer to state variables
     */
    void initializeStates(AmiVector *x);

    /**
     * Initialization of initial state sensitivities
     * @param sx pointer to state variable sensititivies
     * @param x pointer to state variables
     */
    void initializeStateSensitivities(AmiVectorArray *sx, AmiVector *x);

    /**
     * Initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     * @param x pointer to state variables
     * @param dx pointer to time derivative of states (DAE only)
     */
    void initHeaviside(AmiVector *x, AmiVector *dx);

    /**
     * @brief Number of parameters wrt to which sensitivities are computed
     * @return length of sensitivity index vector
     */
    int nplist() const;

    /**
     * @brief Total number of model parameters
     * @return length of parameter vector
     */
    int np() const;

    /**
     * @brief Number of constants
     * @return length of constant vector
     */
    int nk() const;

    /**
     * @brief Number of conservation laws
     * @return difference between nx_rdata and nx_solver
     */
    int ncl() const;

    /**
     * @brief Fixed parameters
     * @return pointer to constants array
     */
    const double *k() const;

    /**
     * @brief Get nmaxevent
     * @return maximum number of events that may occur for each type
     */
    int nMaxEvent() const;

    /**
     * @brief Set nmaxevent
     * @param nmaxevent maximum number of events that may occur for each type
     */
    void setNMaxEvent(int nmaxevent);

    /**
     * @brief Get number of timepoints
     * @return number of timepoints
     */
    int nt() const;

    /**
     * @brief Get ParameterScale for each parameter
     * @return vector of parameter scale
     */
    std::vector<ParameterScaling> const &getParameterScale() const;

    /**
     * @brief Set ParameterScale for each parameter, resets initial state
     * sensitivities
     * @param pscale scalar parameter scale for all parameters
     */
    void setParameterScale(ParameterScaling pscale);

    /**
     * @brief Set ParameterScale for each parameter, resets initial state
     * sensitivities
     * @param pscaleVec vector of parameter scales
     */
    void setParameterScale(const std::vector<ParameterScaling> &pscaleVec);

    /**
     * @brief Get the parameter vector
     * @return The user-set parameters (see also getUnscaledParameters)
     */
    std::vector<realtype> const &getParameters() const;

    /**
     * @brief Sets the parameter vector
     * @param p vector of parameters
     */
    void setParameters(std::vector<realtype> const &p);

    /**
     * @brief Gets parameters with transformation according to ParameterScale
     * applied
     * @return unscaled parameters
     */
    std::vector<realtype> const &getUnscaledParameters() const;

    /**
     * @brief Gets the fixedParameter member
     * @return vector of fixed parameters
     */
    std::vector<realtype> const &getFixedParameters() const;

    /**
     * @brief Sets the fixedParameter member
     * @param k vector of fixed parameters
     */
    void setFixedParameters(std::vector<realtype> const &k);

    /**
     * @brief Get value of fixed parameter with the specified Id
     * @param par_id parameter id
     * @return parameter value
     */
    realtype getFixedParameterById(std::string const &par_id) const;

    /**
     * @brief Get value of fixed parameter with the specified name,
         if multiple parameters have the same name,
         the first parameter with matching name is returned
     * @param par_name parameter name
     * @return parameter value
     */
    realtype getFixedParameterByName(std::string const &par_name) const;

    /**
     * @brief Set value of first fixed parameter with the specified id
     * @param par_id fixed parameter id
     * @param value fixed parameter value
     */
    void setFixedParameterById(std::string const &par_id, realtype value);

    /**
     * @brief Set values of all fixed parameters with the id matching the
     * specified regex
     * @param par_id_regex fixed parameter name regex
     * @param value fixed parameter value
     * @return number of fixed parameter ids that matched the regex
     */
    int setFixedParametersByIdRegex(std::string const &par_id_regex,
                                    realtype value);

    /**
     * @brief Set value of first fixed parameter with the specified name,
     * @param par_name fixed parameter id
     * @param value fixed parameter value
     */
    void setFixedParameterByName(std::string const &par_name, realtype value);

    /**
     * @brief Set value of all fixed parameters with name matching the specified
     * regex,
     * @param par_name_regex fixed parameter name regex
     * @param value fixed parameter value
     * @return number of fixed parameter names that matched the regex
     */
    int setFixedParametersByNameRegex(std::string const &par_name_regex,
                                      realtype value);

    /**
     * @brief Get the timepoint vector
     * @return timepoint vector
     */
    std::vector<realtype> const &getTimepoints() const;

    /**
     * @brief Set the timepoint vector
     * @param ts timepoint vector
     */
    void setTimepoints(std::vector<realtype> const &ts);

    /**
     * @brief gets flags indicating whether states should be treated as
     * non-negative
     * @return vector of flags
     */
    std::vector<bool> const &getStateIsNonNegative() const;

    /**
     * @brief sets flags indicating whether states should be treated as
     * non-negative
     * @param stateIsNonNegative vector of flags
     */
    void setStateIsNonNegative(std::vector<bool> const &stateIsNonNegative);

    /**
     * @brief sets flags indicating that all states should be treated as
     * non-negative
     */
    void setAllStatesNonNegative();

    /**
     * @brief Get timepoint for given index
     * @param idx timepoint index
     * @return timepoint
     */
    double t(int idx) const;

    /**
     * @brief Get the list of parameters for which sensitivities are computed
     * @return list of parameter indices
     */
    std::vector<int> const &getParameterList() const;

    /**
     * @brief Set the list of parameters for which sensitivities are
     * computed, resets initial state sensitivities
     * @param plist list of parameter indices
     */
    void setParameterList(std::vector<int> const &plist);

    /**
     * @brief Get the initial states
     * @return initial state vector
     */
    std::vector<realtype> const &getInitialStates() const;

    /**
     * @brief Set the initial states
     * @param x0 initial state vector
     */
    void setInitialStates(std::vector<realtype> const &x0);

    /**
     * @brief Get the initial states sensitivities
     * @return vector of initial state sensitivities
     */
    std::vector<realtype> const &getInitialStateSensitivities() const;

    /**
     * @brief Set the initial state sensitivities
     * @param sx0 vector of initial state sensitivities with chainrule
     * applied. This could be a slice of ReturnData::sx or ReturnData::sx0
     */
    void setInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief Set the initial state sensitivities
     * @param sx0 vector of initial state sensitivities without chainrule
     * applied. This could be the readin from a model.sx0data saved to hdf5.
     */
    void setUnscaledInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief get simulation start time
     * @return simulation start time
     */
    double t0() const;

    /**
     * @brief set simulation start time
     * @param t0 simulation start time
     */
    void setT0(double t0);

    /**
     * @brief entry in parameter list
     * @param pos index
     * @return entry
     */
    int plist(int pos) const;

    /**
     * @brief Require computation of sensitivities for all parameters p
     * [0..np[ in natural order, resets initial state sensitivities
     */
    void requireSensitivitiesForAllParameters();

    /**
     * @brief Recurring terms in xdot
     * @param t timepoint
     * @param x array with the states
     */
    void fw(const realtype t, const realtype *x);

    /**
     * @brief Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x array with the states
     */
    void fdwdp(const realtype t, const realtype *x);

    /**
     * @brief Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x array with the states
     */
    void fdwdx(const realtype t, const realtype *x);

    /**
     * Residual function
     * @param it time index
     * @param rdata ReturnData instance to which result will be written
     * @param edata ExpData instance containing observable data
     */
    void fres(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Chi-squared function
     * @param it time index
     * @param rdata ReturnData instance to which result will be written
     */
    void fchi2(const int it, ReturnData *rdata);

    /**
     * Residual sensitivity function
     * @param it time index
     * @param rdata ReturnData instance to which result will be written
     * @param edata ExpData instance containing observable data
     */
    void fsres(const int it, ReturnData *rdata, const ExpData *edata);

    /**
     * Fisher information matrix function
     * @param it time index
     * @param rdata ReturnData instance to which result will be written
     */
    void fFIM(const int it, ReturnData *rdata);

    /**
     * Update the heaviside variables h on event occurences
     *
     * @param rootsfound provides the direction of the zero-crossing, so adding
     * it will give the right update to the heaviside variables (zero if no root
     * was found)
     */
    void updateHeaviside(const std::vector<int> &rootsfound);

    /**
     * Updates the heaviside variables h on event occurences in the backward
     problem
     * @param rootsfound provides the direction of the zero-crossing, so adding
         it will give the right update to the heaviside variables (zero if no
     root was found)
     */
    void updateHeavisideB(const int *rootsfound);

    /**
     * @brief Serialize Model (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, Model &r,
                                                const unsigned int version);

    /**
     * @brief Check equality of data members
     * @param a first model instance
     * @param b second model instance
     * @return equality
     */
    friend bool operator==(const Model &a, const Model &b);

    /**
     * Get current timepoint from index
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return current timepoint
     */
    realtype gett(const int it, const ReturnData *rdata) const;

    /**
     * @brief Check if the given array has only finite elements.
     * If not try to give hints by which other fields this could be caused.
     * @param N number of datapoints in array
     * @param array arrays of values
     * @param fun name of the fucntion that generated the values
     * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found,
     * AMICI_SUCCESS otherwise
     */
    int checkFinite(const int N, const realtype *array, const char *fun) const;

    /**
     * @brief Reports whether the model has parameter names set.
     * @return boolean indicating whether parameter names were set
     */
    virtual bool hasParameterNames() const;

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const;

    /**
     * @brief Reports whether the model has state names set.
     * @return boolean indicating whether state names were set
     */
    virtual bool hasStateNames() const;

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const;

    /**
     * @brief Reports whether the model has fixed parameter names set.
     * @return boolean indicating whether fixed parameter names were set
     */
    virtual bool hasFixedParameterNames() const;

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const;

    /**
     * @brief Reports whether the model has observable names set.
     * @return boolean indicating whether observabke names were set
     */
    virtual bool hasObservableNames() const;

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const;

    /**
     * @brief Reports whether the model has parameter ids set.
     * @return boolean indicating whether parameter ids were set
     */
    virtual bool hasParameterIds() const;

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const;

    /**
     * @brief Get value of first model parameter with the specified id
     * @param par_id parameter id
     * @return parameter value
     */
    realtype getParameterById(std::string const &par_id) const;

    /**
     * @brief Get value of first model parameter with the specified name,
     * @param par_name parameter name
     * @return parameter value
     */
    realtype getParameterByName(std::string const &par_name) const;

    /**
     * @brief Set value of first model parameter with the specified id
     * @param par_id parameter id
     * @param value parameter value
     */
    void setParameterById(std::string const &par_id, realtype value);

    /**
     * @brief Set all values of model parameters with ids matching the specified
     * regex
     * @param par_id_regex parameter id regex
     * @param value parameter value
     * @return number of parameter ids that matched the regex
     */
    int setParametersByIdRegex(std::string const &par_id_regex, realtype value);

    /**
     * @brief Set value of first model parameter with the specified name
     * @param par_name parameter name
     * @param value parameter value
     */
    void setParameterByName(std::string const &par_name, realtype value);

    /**
     * @brief Set all values of all model parameters with names matching the
     * specified regex
     * @param par_name_regex parameter name regex
     * @param value parameter value
     * @return number of fixed parameter names that matched the regex
     */
    int setParametersByNameRegex(std::string const &par_name_regex,
                                 realtype value);

    /**
     * @brief Reports whether the model has state ids set.
     * @return
     */
    virtual bool hasStateIds() const;

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const;

    /**
     * @brief Reports whether the model has fixed parameter ids set.
     * @return boolean indicating whether fixed parameter ids were set
     */
    virtual bool hasFixedParameterIds() const;

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getFixedParameterIds() const;

    /**
     * @brief Reports whether the model has observable ids set.
     * @return boolean indicating whether observale ids were set
     */
    virtual bool hasObservableIds() const;

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    virtual std::vector<std::string> getObservableIds() const;

    /**
     * @brief Sets the mode how sensitivities are computed in the steadystate
     * simulation
     * @param mode steadyStateSensitivityMode
     */
    void setSteadyStateSensitivityMode(const SteadyStateSensitivityMode mode);

    /**
     * @brief Gets the mode how sensitivities are computed in the steadystate
     * simulation
     * @return flag value
     */
    SteadyStateSensitivityMode getSteadyStateSensitivityMode() const;

    /**
     * @brief Set whether initial states depending on fixedParmeters are to be
     * reinitialized after preequilibration and presimulation
     * @param flag true/false
     */
    void setReinitializeFixedParameterInitialStates(bool flag);

    /**
     * @brief Get whether initial states depending on fixedParmeters are to be
     * reinitialized after preequilibration and presimulation
     * @return flag true/false
     */
    bool getReinitializeFixedParameterInitialStates() const;

    /**
     * @brief Set if the result of every call to Model::f* should be checked
     * for finiteness
     * @param alwaysCheck
     */
    void setAlwaysCheckFinite(bool alwaysCheck);

    /**
     * @brief Get setting of whether the result of every call to Model::f*
     * should be checked for finiteness
     * @return that
     */
    bool getAlwaysCheckFinite() const;

    /** number of states */
    const int nx_rdata;
    /** number of states in the unaugmented system */
    const int nxtrue_rdata;
    /** number of states with conservation laws applied */
    const int nx_solver;
    /** number of states in the unaugmented system with conservation laws
     * applied */
    const int nxtrue_solver;
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
    const SecondOrderMode o2mode;
    /** index indicating to which event an event output belongs */
    const std::vector<int> z2event;
    /** flag array for DAE equations */
    const std::vector<realtype> idlist;

    /** data standard deviation for current timepoint (dimension: ny) */
    std::vector<realtype> sigmay;
    /** parameter derivative of data standard deviation for current timepoint
     * (dimension: nplist x ny, row-major) */
    std::vector<realtype> dsigmaydp;
    /** event standard deviation for current timepoint (dimension: nz) */
    std::vector<realtype> sigmaz;
    /** parameter derivative of event standard deviation for current timepoint
     * (dimension: nplist x nz, row-major) */
    std::vector<realtype> dsigmazdp;
    /** parameter derivative of data likelihood for current timepoint
     * (dimension: nplist x nJ, row-major) */
    std::vector<realtype> dJydp;
    /** parameter derivative of event likelihood for current timepoint
     * (dimension: nplist x nJ, row-major) */
    std::vector<realtype> dJzdp;

    /** change in x at current timepoint (dimension: nx_solver) */
    std::vector<realtype> deltax;
    /** change in sx at current timepoint (dimension: nplist x nx_solver,
     * row-major) */
    std::vector<realtype> deltasx;
    /** change in xB at current timepoint (dimension: nJ x nxtrue_cl, row-major)
     */
    std::vector<realtype> deltaxB;
    /** change in qB at current timepoint (dimension: nJ x nplist, row-major) */
    std::vector<realtype> deltaqB;

    /** tempory storage of dxdotdp data across functions (dimension: nplist x
     * nx_solver, row-major) */
    std::vector<realtype> dxdotdp;

  protected:
    /**
     * @brief Set the nplist-dependent vectors to their proper sizes
     */
    void initializeVectors();

    /**
     * Model specific implementation of fx_rdata
     * @param x_rdata state variables with conservation laws expanded
     * @param x_solver state variables with conservation laws applied
     * @param tcl total abundances for conservation laws
     **/
    virtual void fx_rdata(realtype *x_rdata, const realtype *x_solver,
                          const realtype *tcl);

    /**
     * Model specific implementation of fsx_solver
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param sx_solver state sensitivity variables with conservation laws
     * applied
     * @param stcl sensitivities of total abundances for conservation laws
     * @param ip sensitivity index
     **/
    virtual void fsx_rdata(realtype *sx_rdata, const realtype *sx_solver,
                           const realtype *stcl, const int ip);

    /**
     * Model specific implementation of fx_solver
     * @param x_solver state variables with conservation laws applied
     * @param x_rdata state variables with conservation laws expanded
     **/
    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata);

    /**
     * Model specific implementation of fsx_solver
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param sx_solver state sensitivity variables with conservation laws
     *applied
     **/
    virtual void fsx_solver(realtype *sx_solver, const realtype *sx_rdata);

    /**
     * Model specific implementation of ftotal_cl
     * @param total_cl total abundances of conservation laws
     * @param x_rdata state variables with conservation laws expanded
     **/
    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata);

    /**
     * Model specific implementation of fstotal_cl
     * @param stotal_cl sensitivites for the total abundances of
     * conservation laws
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param ip sensitivity index
     **/
    virtual void fstotal_cl(realtype *stotal_cl, const realtype *sx_rdata,
                            const int ip);

    /** create my slice at timepoint
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
    void getmy(const int it, const ExpData *edata);

    /** create mz slice at event
     * @param nroots event occurence
     * @param edata pointer to experimental data instance
     */
    void getmz(const int nroots, const ExpData *edata);

    /** create y slice at timepoint
     * @param it timepoint index
     * @param rdata pointer to return data instance
     * @return y y-slice from rdata instance
     */
    const realtype *gety(const int it, const ReturnData *rdata) const;

    /** create z slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return z slice
     */
    const realtype *getz(const int nroots, const ReturnData *rdata) const;

    /** create rz slice at event
     * @param nroots event occurence
     * @param rdata pointer to return data instance
     * @return rz slice
     */
    const realtype *getrz(const int nroots, const ReturnData *rdata) const;

    /** create sz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return z slice
     */
    const realtype *getsz(const int nroots, const int ip,
                          const ReturnData *rdata) const;

    /** create srz slice at event
     * @param nroots event occurence
     * @param ip sensitivity index
     * @param rdata pointer to return data instance
     * @return rz slice
     */
    const realtype *getsrz(const int nroots, const int ip,
                           const ReturnData *rdata) const;

    /**
     * @brief Computes nonnegative state vector according to stateIsNonNegative
     * if anyStateNonNegative is set to false, i.e., all entries in
     * stateIsNonNegative are false, this function directly returns x, otherwise
     * all entries of x are copied in to x_pos_tmp and negative values are
     * replaced by 0 where applicable
     *
     * @param x state vector possibly containing negative values
     * @return state vector with negative values replaced by 0 according to
     * stateIsNonNegative
     */
    N_Vector computeX_pos(N_Vector x);

    /** Sparse Jacobian (dimension: nnz)*/
    SlsMatWrapper J;

    /** current observable (dimension: nytrue) */
    std::vector<realtype> my;
    /** current event measurement (dimension: nztrue) */
    std::vector<realtype> mz;
    /** observable derivative of data likelihood (dimension nJ x nytrue x ny,
     * ordering = ?) */
    std::vector<realtype> dJydy;
    /** observable sigma derivative of data likelihood (dimension nJ x nytrue x
     * ny, ordering = ?) */
    std::vector<realtype> dJydsigma;

    /** event ouput derivative of event likelihood (dimension nJ x nztrue x nz,
     * ordering = ?) */
    std::vector<realtype> dJzdz;
    /** event sigma derivative of event likelihood (dimension nJ x nztrue x nz,
     * ordering = ?) */
    std::vector<realtype> dJzdsigma;
    /** event ouput derivative of event likelihood at final timepoint (dimension
     * nJ x nztrue x nz, ordering = ?) */
    std::vector<realtype> dJrzdz;
    /** event sigma derivative of event likelihood at final timepoint (dimension
     * nJ x nztrue x nz, ordering = ?) */
    std::vector<realtype> dJrzdsigma;
    /** state derivative of event output (dimension: nz x nx_solver, ordering =
     * ?) */
    std::vector<realtype> dzdx;
    /** parameter derivative of event output (dimension: nz x nplist, ordering =
     * ?) */
    std::vector<realtype> dzdp;
    /** state derivative of event timepoint (dimension: nz x nx_solver, ordering
     * = ?) */
    std::vector<realtype> drzdx;
    /** parameter derivative of event timepoint (dimension: nz x nplist,
     * ordering = ?) */
    std::vector<realtype> drzdp;
    /** parameter derivative of observable (dimension: nplist x ny, row-major)
     */
    std::vector<realtype> dydp;

    /** state derivative of observable (dimension: ny x nx_solver, ordering = ?)
     */
    std::vector<realtype> dydx;
    /** tempory storage of w data across functions (dimension: nw) */
    std::vector<realtype> w;
    /** tempory storage of sparse dwdx data across functions (dimension: ndwdx)
     */
    std::vector<realtype> dwdx;
    /** tempory storage of sparse dwdp data across functions (dimension: ndwdp)
     */
    std::vector<realtype> dwdp;
    /** tempory storage of mass matrix data across functions (dimension:
     * nx_solver) */
    std::vector<realtype> M;
    /** tempory storage of stau data across functions (dimension: nplist) */
    std::vector<realtype> stau;

    /** tempory storage of sx data for flattening
         (dimension: nx_solver x nplist, ordering = row-major) */
    std::vector<realtype> sx;

    /** tempory storage x_rdata (dimension: nx_rdata) */
    std::vector<realtype> x_rdata;
    /** tempory storage sx_rdata slice (dimension: nx_rdata) */
    std::vector<realtype> sx_rdata;

    /** flag indicating whether a certain heaviside function should be active or
         not (dimension: ne) */
    std::vector<realtype> h;

    /** unscaled parameters (dimension: np) */
    std::vector<realtype> unscaledParameters;

    /** orignal user-provided, possibly scaled parameter array (size np) */
    std::vector<realtype> originalParameters;

    /** constants (dimension: nk) */
    std::vector<realtype> fixedParameters;

    /** total abundances for conservation laws
         (dimension: nx_rdata-nx_solver) */
    std::vector<realtype> total_cl;

    /** sensitivities of total abundances for conservation laws
         (dimension: (nx_rdata-nx_solver) * np, ordering = row-major)*/
    std::vector<realtype> stotal_cl;

    /** indexes of parameters wrt to which sensitivities are computed (dimension
     * nplist) */
    std::vector<int> plist_;

    /** state initialisation (size nx_solver) */
    std::vector<double> x0data;

    /** sensitivity initialisation (size nx_rdata * nplist, ordering =
     * row-major) */
    std::vector<realtype> sx0data;

    /** timepoints (size nt) */
    std::vector<realtype> ts;

    /** vector of bools indicating whether state variables are to be assumed to
     * be positive */
    std::vector<bool> stateIsNonNegative;

    /** boolean indicating whether any entry in stateIsNonNegative is `true` */
    bool anyStateNonNegative = false;

    /** temporary storage of positified state variables according to
     * stateIsNonNegative */
    AmiVector x_pos_tmp;

    /** maximal number of events to track */
    int nmaxevent = 10;

    /** parameter transformation of `originalParameters` (dimension np) */
    std::vector<ParameterScaling> pscale;

    /** starting time */
    double tstart = 0.0;

    /** flag indicating whether steadystate sensivities are to be computed
     *  via FSA when steadyStateSimulation is used */
    SteadyStateSensitivityMode steadyStateSensitivityMode =
        SteadyStateSensitivityMode::newtonOnly;

    /** flag indicating whether reinitialization of states depending on
     *  fixed parameters is activated
     */
    bool reinitializeFixedParameterInitialStates = false;

    /** Indicates whether the result of every call to Model::f* should be
     * checked for finiteness */
    bool alwaysCheckFinite = false;
};

bool operator==(const Model &a, const Model &b);

} // namespace amici

#endif // AMICI_MODEL_H
