#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

#include "amici/abstract_model.h"
#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <map>
#include <memory>
#include <vector>

namespace amici {

class ExpData;
class Model;
class Solver;
class AmiciApplication;

extern AmiciApplication defaultContext;

} // namespace amici

// for serialization friend in amici::Model
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::Model &u, unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief implements an exchange format to store and transfer the state of the model at a
 * specific timepoint. This is designed to only encompass the minimal number of attributes that need to be
 * transferred.
 */
struct ModelState {
    /** flag indicating whether a certain heaviside function should be active or
         not (dimension: ne) */
    std::vector<realtype> h;

    /** total abundances for conservation laws
     (dimension: nx_rdata-nx_solver) */
    std::vector<realtype> total_cl;

    /** sensitivities of total abundances for conservation laws
     (dimension: (nx_rdata-nx_solver) x np, row-major) */
    std::vector<realtype> stotal_cl;

    /** unscaled parameters (dimension: np) */
    std::vector<realtype> unscaledParameters;

    /** constants (dimension: nk) */
    std::vector<realtype> fixedParameters;

    /** indexes of parameters wrt to which sensitivities are computed
     *  (dimension: nplist)
     */
    std::vector<int> plist;
};



/**
 * @brief The Model class represents an AMICI ODE/DAE model. The model can compute
 * various model related quantities based on symbolically generated code.
 */
class Model : public AbstractModel {
  public:
    /** default constructor */
    Model() = default;

    /**
     * @brief Constructor with model dimensions
     * @param nx_rdata number of state variables
     * @param nxtrue_rdata number of state variables of the non-augmented model
     * @param nx_solver number of state variables with conservation laws applied
     * @param nxtrue_solver number of state variables of the non-augmented model
     * with conservation laws applied
     * @param nx_solver_reinit number of state variables with conservation laws
     * subject to reinitialization
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
     * @param ndxdotdw number of nonzero elements in the w derivative of xdot
     * @param ndJydy number of nonzero elements in the y derivative of dJy
     * (dimension nytrue)
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
    Model(int nx_rdata, int nxtrue_rdata, int nx_solver, int nxtrue_solver,
          int nx_solver_reinit, int ny, int nytrue, int nz, int nztrue, int ne,
          int nJ, int nw, int ndwdx, int ndwdp, int ndxdotdw,
          std::vector<int> ndJydy, int nnz, int ubw, int lbw,
          amici::SecondOrderMode o2mode,
          const std::vector<amici::realtype> &p, std::vector<amici::realtype> k,
          const std::vector<int> &plist, std::vector<amici::realtype> idlist,
          std::vector<int> z2event, bool pythonGenerated = false,
          int ndxdotdp_explicit = 0, int ndxdotdp_implicit = 0);

    /** destructor */
    ~Model() override = default;

    /**
     * @brief Copy assignment is disabled until const members are removed
     * @param other object to copy from
     * @return
     */
    Model &operator=(Model const &other) = delete;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Model *clone() const = 0;

    /**
     * @brief Serialize Model (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param u Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, Model &u,
                                                unsigned int version);

    /**
     * @brief Check equality of data members
     * @param a first model instance
     * @param b second model instance
     * @return equality
     */
    friend bool operator==(const Model &a, const Model &b);

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
    using AbstractModel::fdwdp_colptrs;
    using AbstractModel::fdwdp_rowvals;
    using AbstractModel::fdwdx;
    using AbstractModel::fdwdx_colptrs;
    using AbstractModel::fdwdx_rowvals;
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
     * @brief Initialization of model properties
     * @param x pointer to state variables
     * @param dx pointer to time derivative of states (DAE only)
     * @param sx pointer to state variable sensititivies
     * @param sdx pointer to time derivative of state sensitivities (DAE only)
     * @param computeSensitivities flag indicating whether sensitivities are to
     * be computed
     */
    void initialize(AmiVector &x, AmiVector &dx, AmiVectorArray &sx,
                    AmiVectorArray &sdx, bool computeSensitivities);

    /**
     * @brief Initialization of model properties
     * @param xB adjoint state variables
     * @param dxB time derivative of adjoint states (DAE only)
     * @param xQB adjoint quadratures
     * @param posteq flag indicating whether postequilibration was performed
     */
    void initializeB(AmiVector &xB, AmiVector &dxB, AmiVector &xQB,
                     bool posteq) const;

    /**
     * @brief Initialization of initial states
     * @param x pointer to state variables
     */
    void initializeStates(AmiVector &x);

    /**
     * @brief Initialization of initial state sensitivities
     * @param sx pointer to state variable sensititivies
     * @param x pointer to state variables
     */
    void initializeStateSensitivities(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Initialises the heaviside variables h at the intial time t0
     * heaviside variables activate/deactivate on event occurences
     * @param x pointer to state variables
     * @param dx pointer to time derivative of states (DAE only)
     */
    void initHeaviside(const AmiVector &x, const AmiVector &dx);

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
     * @brief Number of solver states subject to reinitialization
     * @return model member nx_solver_reinit
     */
    int nx_reinit() const;

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
     * @brief Gets parameters with transformation according to ParameterScale
     * applied
     * @return unscaled parameters
     */
    std::vector<realtype> const &getUnscaledParameters() const;

    /**
     * @brief Get the parameter vector
     * @return The user-set parameters (see also getUnscaledParameters)
     */
    std::vector<realtype> const &getParameters() const;

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
     * @brief Sets the parameter vector
     * @param p vector of parameters
     */
    void setParameters(std::vector<realtype> const &p);

    /**
     * @brief Sets model parameters according to the parameter IDs and mapped
     * values.
     * @param p map of parameters IDs and values
     * @param ignoreErrors Ignore errors such as parameter IDs in p which are
     * not model parameters
     */
    void setParameterById(std::map<std::string, realtype> const &p,
                          bool ignoreErrors = false);

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
     * @brief Sets model parameters according to the parameter name and mapped
     * values.
     * @param p map of parameters names and values
     * @param ignoreErrors Ignore errors such as parameter names in p which are
     * not model parameters
     */
    void setParameterByName(std::map<std::string, realtype> const &p,
                            bool ignoreErrors = false);

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
     * @brief Gets the fixedParameter member
     * @return vector of fixed parameters
     */
    std::vector<realtype> const &getFixedParameters() const;

    /**
     * @brief Get value of fixed parameter with the specified Id
     * @param par_id parameter id
     * @return parameter value
     */
    realtype getFixedParameterById(std::string const &par_id) const;

    /**
     * @brief Get value of fixed parameter with the specified name, if multiple
     * parameters have the same name, the first parameter with matching name is
     * returned
     * @param par_name parameter name
     * @return parameter value
     */
    realtype getFixedParameterByName(std::string const &par_name) const;

    /**
     * @brief Sets the fixedParameter member
     * @param k vector of fixed parameters
     */
    void setFixedParameters(std::vector<realtype> const &k);

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
     * @brief Returns the model name
     * @return model name
     */
    virtual std::string getName() const;

    /**
     * @brief Reports whether the model has parameter names set. Also returns
     * true if the number of corresponding variables is just zero.
     * @return boolean indicating whether parameter names were set
     */
    virtual bool hasParameterNames() const;

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const;

    /**
     * @brief Reports whether the model has state names set. Also returns true
     * if the number of corresponding variables is just zero.
     * @return boolean indicating whether state names were set
     */
    virtual bool hasStateNames() const;

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const;

    /**
     * @brief Reports whether the model has fixed parameter names set. Also
     * returns true if the number of corresponding variables is just zero.
     * @return boolean indicating whether fixed parameter names were set
     */
    virtual bool hasFixedParameterNames() const;

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const;

    /**
     * @brief Reports whether the model has observable names set. Also returns
     * true if the number of corresponding variables is just zero.
     * @return boolean indicating whether observabke names were set
     */
    virtual bool hasObservableNames() const;

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const;

    /**
     * @brief Reports whether the model has parameter ids set. Also returns true
     * if the number of corresponding variables is just zero.
     * @return boolean indicating whether parameter ids were set
     */
    virtual bool hasParameterIds() const;

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getParameterIds() const;

    /**
     * @brief Reports whether the model has state ids set. Also returns true if
     * the number of corresponding variables is just zero.
     * @return boolean indicating whether state ids were set
     */
    virtual bool hasStateIds() const;

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    virtual std::vector<std::string> getStateIds() const;

    /**
     * @brief Reports whether the model has fixed parameter ids set. Also
     * returns true if the number of corresponding variables is just zero.
     * @return boolean indicating whether fixed parameter ids were set
     */
    virtual bool hasFixedParameterIds() const;

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    virtual std::vector<std::string> getFixedParameterIds() const;

    /**
     * @brief Reports whether the model has observable ids set. Also returns
     * true if the number of corresponding variables is just zero.
     * @return boolean indicating whether observale ids were set
     */
    virtual bool hasObservableIds() const;

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    virtual std::vector<std::string> getObservableIds() const;

    /**
     * @brief Get the timepoint vector
     * @return timepoint vector
     */
    std::vector<realtype> const &getTimepoints() const;

    /**
     * @brief get simulation timepoint for time index it
     * @param it time index
     * @return t timepoint
     */
    realtype getTimepoint(int it) const;

    /**
     * @brief Set the timepoint vector
     * @param ts timepoint vector
     */
    void setTimepoints(std::vector<realtype> const &ts);

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
     * @brief retrieves the current model state
     * @return current model state
     */
    ModelState const &getModelState() const {
        return state_;
    };

    /**
     * @brief sets the current model state
     * @param state model state
     */
    void setModelState(ModelState const &state) {
        if (static_cast<int>(state.unscaledParameters.size()) != np())
            throw AmiException("Mismatch in parameter size");
        if (static_cast<int>(state.fixedParameters.size()) != nk())
            throw AmiException("Mismatch in fixed parameter size");
        if (static_cast<int>(state.h.size()) != ne)
            throw AmiException("Mismatch in heaviside size");
        if (static_cast<int>(state.total_cl.size()) != ncl())
            throw AmiException("Mismatch in conservation law size");
        if (static_cast<int>(state.stotal_cl.size()) != ncl() * np() )
            throw AmiException("Mismatch in conservation law sensitivity size");
        state_ = state;
    };

    /**
     * @brief Get the list of parameters for which sensitivities are computed
     * @return list of parameter indices
     */
    std::vector<int> const &getParameterList() const;

    /**
     * @brief entry in parameter list
     * @param pos index
     * @return entry
     */
    int plist(int pos) const;

    /**
     * @brief Set the list of parameters for which sensitivities are computed,
     * resets initial state sensitivities
     * @param plist list of parameter indices
     */
    void setParameterList(std::vector<int> const &plist);

    /**
     * @brief Get the initial states
     * @return initial state vector
     */
    std::vector<realtype> getInitialStates();

    /**
     * @brief Set the initial states
     * @param x0 initial state vector
     */
    void setInitialStates(std::vector<realtype> const &x0);

    /**
     * @brief Returns whether custom initial states have been set
     * @return True if has custom initial states, otherwise false
     */
    bool hasCustomInitialStates() const;

    /**
     * @brief Get the initial states sensitivities
     * @return vector of initial state sensitivities
     */
    std::vector<realtype> getInitialStateSensitivities();

    /**
     * @brief Set the initial state sensitivities
     * @param sx0 vector of initial state sensitivities with chainrule applied.
     * This could be a slice of ReturnData::sx or ReturnData::sx0
     */
    void setInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief Returns whether custom initial state sensitivitiess have been set
     * @return True if has custom initial state sensitivitiess, otherwise false
     */
    bool hasCustomInitialStateSensitivities() const;

    /**
     * @brief Set the initial state sensitivities
     * @param sx0 vector of initial state sensitivities without chainrule
     * applied. This could be the readin from a model.sx0data saved to hdf5.
     */
    void setUnscaledInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief Sets the mode how sensitivities are computed in the steadystate
     * simulation
     * @param mode steadyStateSensitivityMode
     */
    void setSteadyStateSensitivityMode(SteadyStateSensitivityMode mode);

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
     * @brief Require computation of sensitivities for all parameters p [0..np[
     * in natural order, resets initial state sensitivities
     */
    void requireSensitivitiesForAllParameters();

    /**
     * @brief Time-resolved w,
     * @param w buffer (dimension: nw)
     * @param t current timepoint
     * @param x current state
     */
    void getExpression(gsl::span<realtype> w, const realtype t, const AmiVector &x);

    /**
     * @brief Time-resolved observables,
     * @param y buffer (dimension: ny)
     * @param t current timepoint
     * @param x current state
     */
    void getObservable(gsl::span<realtype> y, const realtype t,
                       const AmiVector &x);

    /**
     * @brief Sensitivity of time-resolved observables, total derivative sy =
     * dydx * sx + dydp (only for forward sensitivities)
     * @param sy buffer (dimension: ny x nplist, row-major)
     * @param t timpoint
     * @param x state variables
     * @param sx state sensitivities
     */
    void getObservableSensitivity(gsl::span<realtype> sy, const realtype t,
                                  const AmiVector &x, const AmiVectorArray &sx);

    /**
     * @brief Time-resolved observable standard deviations
     * @param sigmay buffer (dimension: ny)
     * @param it timepoint index
     * @param edata pointer to experimental data instance (optional, pass
     * nullptr to ignore)
     */
    void getObservableSigma(gsl::span<realtype> sigmay, const int it,
                            const ExpData *edata);

    /**
     * @brief Sensitivity of time-resolved observable standard deviation, total
     * derivative (can be used with both adjoint and forward sensitivity)
     * @param ssigmay buffer (dimension: ny x nplist, row-major)
     * @param it timepoint index
     * @param edata pointer to experimental data instance (optional, pass
     * nullptr to ignore)
     */
    void getObservableSigmaSensitivity(gsl::span<realtype> ssigmay,
                                       const int it, const ExpData *edata);

    /**
     * @brief Time-resolved measurement negative log-likelihood Jy
     * @param Jy buffer (dimension: 1)
     * @param it timepoint index
     * @param x state variables
     * @param edata experimental data instance
     */
    void addObservableObjective(realtype &Jy, const int it, const AmiVector &x,
                                const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy, total derivative (to be used with forward senstivities)
     * @param sllh first order buffer (dimension: nplist)
     * @param s2llh second order buffer (dimension: nJ-1 x nplist, row-major)
     * @param it timepoint index
     * @param x state variables
     * @param sx state sensitivities
     * @param edata experimental data instance
     */
    void addObservableObjectiveSensitivity(std::vector<realtype> &sllh,
                                           std::vector<realtype> &s2llh,
                                           const int it, const AmiVector &x,
                                           const AmiVectorArray &sx,
                                           const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy, partial derivative (to be used with adjoint senstivities)
     * @param sllh first order buffer (dimension: nplist)
     * @param s2llh second order buffer (dimension: nJ-1 x nplist, row-major)
     * @param it timepoint index
     * @param x state variables
     * @param edata experimental data instance
     */
    void addPartialObservableObjectiveSensitivity(std::vector<realtype> &sllh,
                                                  std::vector<realtype> &s2llh,
                                                  const int it,
                                                  const AmiVector &x,
                                                  const ExpData &edata);

    /**
     * @brief State sensitivity of the negative loglikelihood Jy, partial
     * derivative (to be used with adjoint senstivities)
     * @param dJydx buffer (dimension: nJ x nx_solver, row-major)
     * @param it timepoint index
     * @param x state variables
     * @param edata experimental data instance
     */
    void getAdjointStateObservableUpdate(gsl::span<realtype> dJydx,
                                         const int it, const AmiVector &x,
                                         const ExpData &edata);

    /**
     * @brief Event-resolved observables
     * @param z buffer (dimension: nz)
     * @param ie event index
     * @param t timepoint
     * @param x state variables
     */
    void getEvent(gsl::span<realtype> z, const int ie, const realtype t,
                  const AmiVector &x);
    /**
     * @brief Sensitivities of event-resolved observables, total derivative,
     * total derivative (only forward sensitivities)
     * @param sz buffer (dimension: nz x nplist, row-major)
     * @param ie event index
     * @param t timepoint
     * @param x state variables
     * @param sx state sensitivities
     */
    void getEventSensitivity(gsl::span<realtype> sz, const int ie,
                             const realtype t, const AmiVector &x,
                             const AmiVectorArray &sx);

    /**
     * @brief Sensitivity of z at final timepoint (ignores sensitivity of
     * timepoint), total derivative
     * @param sz output buffer (dimension: nz x nplist, row-major)
     * @param ie event index
     */
    void getUnobservedEventSensitivity(gsl::span<realtype> sz, const int ie);

    /**
     * @brief Regularization for event-resolved observables
     * @param rz buffer (dimension: nz)
     * @param ie event index
     * @param t timepoint
     * @param x state variables
     */
    void getEventRegularization(gsl::span<realtype> rz, const int ie,
                                const realtype t, const AmiVector &x);

    /**
     * @brief Sensitivities of regularization for event-resolved observables,
     * total derivative (only forward sensitivities)
     * @param srz buffer (dimension: nz x nplist, row-major)
     * @param ie event index
     * @param t timepoint
     * @param x state variables
     * @param sx state sensitivities
     */
    void getEventRegularizationSensitivity(gsl::span<realtype> srz,
                                           const int ie, const realtype t,
                                           const AmiVector &x,
                                           const AmiVectorArray &sx);
    /**
     * @brief Event-resolved observable standard deviations
     * @param sigmaz buffer (dimension: nz)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param edata pointer to experimental data instance (optional, pass
     * nullptr to ignore)
     */
    void getEventSigma(gsl::span<realtype> sigmaz, const int ie,
                       const int nroots, const realtype t,
                       const ExpData *edata);

    /**
     * @brief Sensitivities of event-resolved observable standard deviations,
     * total derivative (only forward sensitivities)
     * @param ssigmaz buffer (dimension: nz x nplist, row-major)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param edata pointer to experimental data instance (optional, pass
     * nullptr to ignore)
     */
    void getEventSigmaSensitivity(gsl::span<realtype> ssigmaz, const int ie,
                                  const int nroots, const realtype t,
                                  const ExpData *edata);

    /**
     * @brief Event-resolved observable negative log-likelihood,
     * @param Jz buffer (dimension: 1)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param x state variables
     * @param edata experimental data instance
     */
    void addEventObjective(realtype &Jz, const int ie, const int nroots,
                           const realtype t, const AmiVector &x,
                           const ExpData &edata);

    /**
     * @brief Event-resolved observable negative log-likelihood,
     * @param Jrz buffer (dimension: 1)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param x state variables
     * @param edata experimental data instance
     */
    void addEventObjectiveRegularization(realtype &Jrz, const int ie,
                                         const int nroots, const realtype t,
                                         const AmiVector &x,
                                         const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy, total derivative (to be used with forward senstivities)
     * @param sllh first order buffer (dimension: nplist)
     * @param s2llh second order buffer (dimension: nJ-1 x nplist, row-major)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param x state variables
     * @param sx state sensitivities
     * @param edata experimental data instance
     */
    void addEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                      std::vector<realtype> &s2llh,
                                      const int ie, const int nroots,
                                      const realtype t, const AmiVector &x,
                                      const AmiVectorArray &sx,
                                      const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy, partial derivative (to be used with adjoint senstivities)
     * @param sllh first order buffer (dimension: nplist)
     * @param s2llh second order buffer (dimension: nJ-1 x nplist, row-major)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param x state variables
     * @param edata experimental data instance
     */
    void addPartialEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                             std::vector<realtype> &s2llh,
                                             const int ie, const int nroots,
                                             const realtype t,
                                             const AmiVector &x,
                                             const ExpData &edata);

    /**
     * @brief State sensitivity of the negative loglikelihood Jz, partial
     * derivative (to be used with adjoint senstivities)
     * @param dJzdx buffer (dimension: nJ x nx_solver, row-major)
     * @param ie event index
     * @param nroots event occurence
     * @param t timepoint
     * @param x state variables
     * @param edata experimental data instance
     */
    void getAdjointStateEventUpdate(gsl::span<realtype> dJzdx, const int ie,
                                    const int nroots, const realtype t,
                                    const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of event timepoint, total derivative (only forward
     * sensi)
     * @param stau current timepoint (dimension: nplist)
     * @param t timepoint
     * @param ie event index
     * @param x state variables
     * @param sx state sensitivities
     */
    void getEventTimeSensitivity(std::vector<realtype> &stau, const realtype t,
                                 const int ie, const AmiVector &x,
                                 const AmiVectorArray &sx);

    /**
     * @brief Update state variables after event
     * @param x current state (will be overwritten)
     * @param ie event index
     * @param t current timepoint
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void addStateEventUpdate(AmiVector &x, const int ie, const realtype t,
                             const AmiVector &xdot, const AmiVector &xdot_old);

    /**
     * @brief Update state sensitivity after event
     * @param sx current state sensitivity (will be overwritten)
     * @param ie event index
     * @param t current timepoint
     * @param x_old current state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     * @param stau timepoint sensitivity, to be computed with
     * Model::getEventTimeSensitivity
     */
    void addStateSensitivityEventUpdate(AmiVectorArray &sx, const int ie,
                                        const realtype t,
                                        const AmiVector &x_old,
                                        const AmiVector &xdot,
                                        const AmiVector &xdot_old,
                                        const std::vector<realtype> &stau);

    /**
     * @brief Update adjoint state after event
     * @param xB current adjoint state (will be overwritten)
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void addAdjointStateEventUpdate(AmiVector &xB, const int ie,
                                    const realtype t, const AmiVector &x,
                                    const AmiVector &xdot,
                                    const AmiVector &xdot_old);

    /**
     * @brief Update adjoint quadratures after event
     * @param xQB current quadrature state (will be overwritten)
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     * @param xB current adjoint state
     * @param xdot current residual function values
     * @param xdot_old value of residual function before event
     */
    void addAdjointQuadratureEventUpdate(AmiVector xQB, const int ie,
                                         const realtype t, const AmiVector &x,
                                         const AmiVector &xB,
                                         const AmiVector &xdot,
                                         const AmiVector &xdot_old);

    /**
     * @brief Update the heaviside variables h on event occurences
     *
     * @param rootsfound provides the direction of the zero-crossing, so adding
     * it will give the right update to the heaviside variables (zero if no root
     * was found)
     */
    void updateHeaviside(const std::vector<int> &rootsfound);

    /**
     * @brief Updates the heaviside variables h on event occurences in the
     * backward problem
     * @param rootsfound provides the direction of the zero-crossing, so adding
     * it will give the right update to the heaviside variables (zero if no root
     * was found)
     */
    void updateHeavisideB(const int *rootsfound);

    /**
     * @brief Check if the given array has only finite elements. If not try to
     * give hints by which other fields this could be caused.
     *
     * @param array arrays of values
     * @param fun name of the fucntion that generated the values
     * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found,
     * AMICI_SUCCESS otherwise
     */
    int checkFinite(gsl::span<const realtype> array, const char *fun) const;

    /**
     * @brief Set if the result of every call to Model::f* should be checked for
     * finiteness
     * @param alwaysCheck
     */
    void setAlwaysCheckFinite(bool alwaysCheck);

    /**
     * @brief Get setting of whether the result of every call to Model::f*
     * should be checked for finiteness
     * @return that
     */
    bool getAlwaysCheckFinite() const;

    /**
     * @brief Initial states
     * @param x pointer to state variables
     */
    void fx0(AmiVector &x);

    /**
     * @brief Sets only those initial states that are specified via
     * fixedParmeters
     * @param x pointer to state variables
     */
    void fx0_fixedParameters(AmiVector &x);

    /**
     * @brief Initial value for initial state sensitivities
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     */
    void fsx0(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Sets only those initial states sensitivities that are affected
     * from fx0 fixedParmeters
     * @param sx pointer to state sensitivity variables
     * @param x pointer to state variables
     */
    void fsx0_fixedParameters(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Sensitivity of derivative initial states sensitivities sdx0 (only
     *  necessary for DAEs)
     */
    virtual void fsdx0();

    /**
     * @brief Expands conservation law for states
     * @param x_rdata pointer to state variables with conservation laws expanded
     * (stored in rdata)
     * @param x_solver pointer to state variables with conservation laws applied
     * (solver returns this)
     */
    void fx_rdata(AmiVector &x_rdata, const AmiVector &x_solver);

    /**
     * @brief Expands conservation law for state sensitivities
     * @param sx_rdata pointer to state variable sensitivities with conservation
     * laws expanded (stored in rdata)
     * @param sx_solver pointer to state variable sensitivities with
     * conservation laws applied (solver returns this)
     */
    void fsx_rdata(AmiVectorArray &sx_rdata, const AmiVectorArray &sx_solver);

    /** number of states */
    int nx_rdata{0};

    /** number of states in the unaugmented system */
    int nxtrue_rdata{0};

    /** number of states with conservation laws applied */
    int nx_solver{0};

    /**
     * number of states in the unaugmented system with conservation laws
     * applied
     */
    int nxtrue_solver{0};

    /** number of solver states subject to reinitialization */
    int nx_solver_reinit{0};

    /** number of observables */
    int ny{0};

    /** number of observables in the unaugmented system */
    int nytrue{0};

    /** number of event outputs */
    int nz{0};

    /** number of event outputs in the unaugmented system */
    int nztrue{0};

    /** number of events */
    int ne{0};

    /** number of common expressions */
    int nw{0};

    /** number of derivatives of common expressions wrt x */
    int ndwdx{0};

    /** number of derivatives of common expressions wrt p */
    int ndwdp{0};

    /** number of nonzero entries in dxdotdw */
    int ndxdotdw{0};
    /** number of nonzero entries in dJydy */
    std::vector<int> ndJydy;

    /** number of nonzero entries in jacobian */
    int nnz{0};

    /** dimension of the augmented objective function for 2nd order ASA */
    int nJ{0};

    /** upper bandwith of the jacobian */
    int ubw{0};

    /** lower bandwith of the jacobian */
    int lbw{0};

    /** flag indicating Matlab or python based model generation */
    bool pythonGenerated;

    /** number of nonzero entries in ndxdotdp_explicit */
    int ndxdotdp_explicit{0};

    /** number of nonzero entries in ndxdotdp_implicit */
    int ndxdotdp_implicit{0};

    /**
     * flag indicating whether for sensi == AMICI_SENSI_ORDER_SECOND
     * directional or full second order derivative will be computed
     */
    SecondOrderMode o2mode{SecondOrderMode::none};

    /** flag array for DAE equations */
    std::vector<realtype> idlist;

    /**
     * temporary storage of dxdotdp data across functions, Python only
     * (dimension: nplist x nx_solver, nnz: ndxdotdp_explicit, type CSC_MAT)
     */
    mutable SUNMatrixWrapper dxdotdp_explicit;

    /**
     * temporary storage of dxdotdp_implicit data across functions, Python only
     * (dimension: nplist x * nx_solver, nnz: ndxdotdp_implicit, type CSC_MAT)
     */
    mutable SUNMatrixWrapper dxdotdp_implicit;

    /**
     * temporary storage of dxdotdp data across functions, Matlab only
     * (dimension: nplist x nx_solver, row-major)
     */
    AmiVectorArray dxdotdp {0, 0};
    /** AMICI context */
    AmiciApplication *app = &defaultContext;

  protected:
    /**
     * @brief Writes part of a slice to a buffer according to indices specified
     * in z2event
     * @param slice source data slice
     * @param buffer output data slice
     * @param ie event index
     */
    void writeSliceEvent(gsl::span<const realtype> slice,
                         gsl::span<realtype> buffer, const int ie);

    /**
     * @brief Writes part of a sensitivity slice to a buffer according to
     * indices specified in z2event
     * @param slice source data slice
     * @param buffer output data slice
     * @param ie event index
     */
    void writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                    gsl::span<realtype> buffer, const int ie);

    /**
     * @brief Seperates first and second order objective sensitivity information
     * and writes them into the respective buffers
     * @param dLLhdp data with mangled first and second order information
     * @param sllh first order buffer
     * @param s2llh second order buffer
     */
    void writeLLHSensitivitySlice(const std::vector<realtype> &dLLhdp,
                                  std::vector<realtype> &sllh,
                                  std::vector<realtype> &s2llh);

    /**
     * @brief Verifies that the provided buffers have the expected size.
     * @param sllh first order buffer
     * @param s2llh second order buffer
     */
    void checkLLHBufferSize(const std::vector<realtype> &sllh,
                            const std::vector<realtype> &s2llh) const;

    /**
     * @brief Set the nplist-dependent vectors to their proper sizes
     */
    void initializeVectors();

    /**
     * @brief Observables / measurements
     * @param t current timepoint
     * @param x current state
     */
    void fy(realtype t, const AmiVector &x);

    /**
     * @brief Partial derivative of observables y w.r.t. model parameters p
     * @param t current timepoint
     * @param x current state
     */
    void fdydp(realtype t, const AmiVector &x);

    /**
     * @brief Partial derivative of observables y w.r.t. state variables x
     * @param t current timepoint
     * @param x current state
     */
    void fdydx(realtype t, const AmiVector &x);

    /**
     * @brief Standard deviation of measurements
     * @param it timepoint index
     * @param edata pointer to experimental data instance
     */
    void fsigmay(int it, const ExpData *edata);

    /**
     * @brief Partial derivative of standard deviation of measurements w.r.t.
     * model
     * @param it timepoint index
     * @param edata pointer to ExpData data instance holding sigma values
     */
    void fdsigmaydp(int it, const ExpData *edata);

    /**
     * @brief Negative log-likelihood of measurements y
     * @param Jy variable to which llh will be added
     * @param it timepoint index
     * @param y simulated observable
     * @param edata pointer to experimental data instance
     */
    void fJy(realtype &Jy, int it, const AmiVector &y, const ExpData &edata);

    /**
     * @brief Model specific implementation of fdJydy colptrs
     * @param indexptrs column pointers
     * @param index ytrue index
     */
    virtual void fdJydy_colptrs(sunindextype *indexptrs, int index);

    /**
     * @brief Model specific implementation of fdJydy row vals
     * @param indexptrs row val pointers
     * @param index ytrue index
     */
    virtual void fdJydy_rowvals(sunindextype *indexptrs, int index);

    /**
     * @brief Partial derivative of time-resolved measurement negative
     * log-likelihood Jy
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydy(int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy w.r.t. standard deviation sigma
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydsigma(int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute sensitivity of time-resolved measurement negative
     * log-likelihood Jy w.r.t. parameters for the given timepoint. Add result
     * to respective fields in rdata.
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydp(const int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy w.r.t. state variables
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydx(const int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Event-resolved output
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fz(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Partial derivative of event-resolved output z w.r.t. to model
     * parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdzdp(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Partial derivative of event-resolved output z w.r.t. to model
     * states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdzdx(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Event root function of events (equal to froot but does not include
     * non-output events)
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void frz(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Sensitivity of event-resolved root output w.r.t. to model
     * parameters p
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdrzdp(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Sensitivity of event-resolved measurements rz w.r.t. to model
     * states x
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdrzdx(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Standard deviation of events
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param edata pointer to experimental data instance
     */
    void fsigmaz(const int ie, const int nroots, const realtype t,
                 const ExpData *edata);

    /**
     * @brief Sensitivity of standard deviation of events measurements w.r.t.
     * model parameters p
     * @param ie event index
     * @param nroots event occurence
     * @param t current timepoint
     * @param edata pointer to experimental data instance
     */
    void fdsigmazdp(int ie, int nroots, realtype t, const ExpData *edata);

    /**
     * @brief Negative log-likelihood of event-resolved measurements z
     * @param Jz variable to which llh will be added
     * @param nroots event index
     * @param z simulated event
     * @param edata pointer to experimental data instance
     */
    void fJz(realtype &Jz, int nroots, const AmiVector &z,
             const ExpData &edata);

    /**
     * @brief Partial derivative of event measurement negative log-likelihood Jz
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJzdz(const int ie, const int nroots, const realtype t,
                const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of event measurement negative log-likelihood Jz w.r.t.
     * standard deviation sigmaz
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJzdsigma(const int ie, const int nroots, const realtype t,
                    const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of event-resolved measurement negative log-likelihood
     * Jz w.r.t. parameters
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJzdp(const int ie, const int nroots, realtype t, const AmiVector &x,
                const ExpData &edata);

    /**
     * @brief Sensitivity of event-resolved measurement negative log-likelihood
     * Jz w.r.t. state variables
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJzdx(const int ie, const int nroots, realtype t, const AmiVector &x,
                const ExpData &edata);

    /**
     * @brief Regularization of negative log-likelihood with roots of
     * event-resolved measurements rz
     * @param Jrz variable to which regularization will be added
     * @param nroots event index
     * @param rz regularization variable
     * @param edata pointer to experimental data instance
     */
    void fJrz(realtype &Jrz, int nroots, const AmiVector &rz,
              const ExpData &edata);

    /**
     * @brief Partial derivative of event measurement negative log-likelihood Jz
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJrzdz(const int ie, const int nroots, const realtype t,
                 const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of event measurement negative log-likelihood Jz w.r.t.
     * standard deviation sigmaz
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJrzdsigma(const int ie, const int nroots, const realtype t,
                     const AmiVector &x, const ExpData &edata);

    /**
     * @brief Recurring terms in xdot
     * @param t timepoint
     * @param x array with the states
     */
    void fw(realtype t, const realtype *x);

    /**
     * @brief Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x array with the states
     */
    void fdwdp(realtype t, const realtype *x);

    /**
     * @brief Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x array with the states
     */
    void fdwdx(realtype t, const realtype *x);

    /**
     * @brief Model specific implementation of fx_rdata
     * @param x_rdata state variables with conservation laws expanded
     * @param x_solver state variables with conservation laws applied
     * @param tcl total abundances for conservation laws
     */
    virtual void fx_rdata(realtype *x_rdata, const realtype *x_solver,
                          const realtype *tcl);

    /**
     * @brief Model specific implementation of fsx_solver
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param sx_solver state sensitivity variables with conservation laws
     * applied
     * @param stcl sensitivities of total abundances for conservation laws
     * @param ip sensitivity index
     */
    virtual void fsx_rdata(realtype *sx_rdata, const realtype *sx_solver,
                           const realtype *stcl, int ip);

    /**
     * @brief Model specific implementation of fx_solver
     * @param x_solver state variables with conservation laws applied
     * @param x_rdata state variables with conservation laws expanded
     */
    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata);

    /**
     * @brief Model specific implementation of fsx_solver
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param sx_solver state sensitivity variables with conservation laws
     * applied
     */
    virtual void fsx_solver(realtype *sx_solver, const realtype *sx_rdata);

    /**
     * @brief Model specific implementation of ftotal_cl
     * @param total_cl total abundances of conservation laws
     * @param x_rdata state variables with conservation laws expanded
     */
    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata);

    /**
     * @brief Model specific implementation of fstotal_cl
     * @param stotal_cl sensitivites for the total abundances of conservation
     * laws
     * @param sx_rdata state sensitivity variables with conservation laws
     * expanded
     * @param ip sensitivity index
     */
    virtual void fstotal_cl(realtype *stotal_cl, const realtype *sx_rdata,
                            int ip);

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
    N_Vector computeX_pos(const_N_Vector x);

    /** all variables necessary for function evaluation */
    ModelState state_;

    /** Sparse Jacobian (dimension: nnz)*/
    mutable SUNMatrixWrapper J_;

    /** Sparse dxdotdw temporary storage (dimension: ndxdotdw) */
    mutable SUNMatrixWrapper dxdotdw_;

    /** Sparse dwdp temporary storage (dimension: ndwdp) */
    mutable SUNMatrixWrapper dwdp_;

    /** Sparse dwdx temporary storage (dimension: ndwdx) */
    mutable SUNMatrixWrapper dwdx_;

    /** Dense Mass matrix (dimension: nx_solver x nx_solver) */
    mutable SUNMatrixWrapper M_;

    /** current observable (dimension: nytrue) */
    mutable std::vector<realtype> my_;

    /** current event measurement (dimension: nztrue) */
    mutable std::vector<realtype> mz_;

    /** Sparse observable derivative of data likelihood, only used if
     * pythonGenerated==true (dimension nytrue, nJ x ny, row-major) */
    mutable std::vector<SUNMatrixWrapper> dJydy_;

    /** observable derivative of data likelihood, only used if
     * pythonGenerated==false (dimension nJ x ny x nytrue, row-major)
     */
    mutable std::vector<realtype> dJydy_matlab_;

    /** observable sigma derivative of data likelihood
     * (dimension nJ x ny x nytrue, row-major)
     */
    mutable std::vector<realtype> dJydsigma_;

    /** state derivative of data likelihood
     * (dimension nJ x nx_solver, row-major)
     */
    mutable std::vector<realtype> dJydx_;

    /** parameter derivative of data likelihood for current timepoint
     * (dimension: nJ x nplist, row-major)
     */
    mutable std::vector<realtype> dJydp_;

    /** event output derivative of event likelihood
     * (dimension nJ x nz x nztrue, row-major)
     */
    mutable std::vector<realtype> dJzdz_;

    /** event sigma derivative of event likelihood
     * (dimension nJ x nz x nztrue, row-major)
     */
    mutable std::vector<realtype> dJzdsigma_;

    /** event output derivative of event likelihood at final timepoint
     * (dimension nJ x nz x nztrue, row-major)
     */
    mutable std::vector<realtype> dJrzdz_;

    /** event sigma derivative of event likelihood at final timepoint
     * (dimension nJ x nz x nztrue, row-major)
     */
    mutable std::vector<realtype> dJrzdsigma_;

    /** state derivative of event likelihood
     * (dimension nJ x nx_solver, row-major)
     */
    mutable std::vector<realtype> dJzdx_;

    /** parameter derivative of event likelihood for current timepoint
     * (dimension: nJ x nplist x, row-major)
     */
    mutable std::vector<realtype> dJzdp_;

    /** state derivative of event output
     * (dimension: nz x nx_solver, row-major)
     */
    mutable std::vector<realtype> dzdx_;

    /** parameter derivative of event output
     * (dimension: nz x nplist, row-major)
     */
    mutable std::vector<realtype> dzdp_;

    /** state derivative of event regularization variable
     * (dimension: nz x nx_solver, row-major)
     */
    mutable std::vector<realtype> drzdx_;

    /** parameter derivative of event regularization variable
     * (dimension: nz x nplist, row-major)
     */
    mutable std::vector<realtype> drzdp_;

    /** parameter derivative of observable
     * (dimension: ny x nplist, row-major)
     */
    mutable std::vector<realtype> dydp_;

    /** state derivative of time-resolved observable
     * (dimension: nx_solver x ny, row-major)
     */
    mutable std::vector<realtype> dydx_;

    /** tempory storage of w data across functions (dimension: nw) */
    mutable std::vector<realtype> w_;

    /** tempory storage for flattened sx,
     * (dimension: nx_solver x nplist, row-major)
     */
    mutable std::vector<realtype> sx_;

    /** tempory storage for x_rdata (dimension: nx_rdata) */
    mutable std::vector<realtype> x_rdata_;

    /** tempory storage for sx_rdata slice (dimension: nx_rdata) */
    mutable std::vector<realtype> sx_rdata_;

    /** temporary storage for time-resolved observable (dimension: ny) */
    mutable std::vector<realtype> y_;

    /** data standard deviation for current timepoint (dimension: ny) */
    mutable std::vector<realtype> sigmay_;

    /** temporary storage for  parameter derivative of data standard deviation,
     * (dimension: ny x nplist, row-major)
     */
    mutable std::vector<realtype> dsigmaydp_;

    /** temporary storage for event-resolved observable (dimension: nz) */
    mutable std::vector<realtype> z_;

    /** temporary storage for event regularization (dimension: nz) */
    mutable std::vector<realtype> rz_;

    /** temporary storage for event standard deviation (dimension: nz) */
    mutable std::vector<realtype> sigmaz_;

    /** temporary storage for  parameter derivative of event standard deviation,
     * (dimension: nz x nplist, row-major)
     */
    mutable std::vector<realtype> dsigmazdp_;

    /** temporary storage for change in x after event (dimension: nx_solver) */
    mutable std::vector<realtype> deltax_;

    /** temporary storage for change in sx after event
     * (dimension: nx_solver x nplist, row-major)
     */
    mutable std::vector<realtype> deltasx_;

    /** temporary storage for change in xB after event (dimension: nx_solver) */
    mutable std::vector<realtype> deltaxB_;

    /** temporary storage for change in qB after event
     * (dimension: nJ x nplist, row-major)
     */
    mutable std::vector<realtype> deltaqB_;

    /** temporary storage of positified state variables according to
     * stateIsNonNegative (dimension: nx_solver) */
    mutable AmiVector x_pos_tmp_ {0};

    /** orignal user-provided, possibly scaled parameter array (dimension: np)
     */
    std::vector<realtype> original_parameters_;

    /** index indicating to which event an event output belongs */
    std::vector<int> z2event_;

    /** state initialisation (size nx_solver) */
    std::vector<realtype> x0data_;

    /** sensitivity initialisation (size nx_rdata x nplist, row-major) */
    std::vector<realtype> sx0data_;

    /** timepoints (size nt) */
    std::vector<realtype> ts_;

    /** vector of bools indicating whether state variables are to be assumed to
     * be positive */
    std::vector<bool> state_is_non_negative_;

    /** boolean indicating whether any entry in stateIsNonNegative is `true` */
    bool any_state_non_negative_ {false};

    /** maximal number of events to track */
    int nmaxevent_ {10};

    /** parameter transformation of `originalParameters` (dimension np) */
    std::vector<ParameterScaling> pscale_;

    /** starting time */
    realtype tstart_ {0.0};

    /** flag indicating whether steadystate sensivities are to be computed
     *  via FSA when steadyStateSimulation is used
     */
    SteadyStateSensitivityMode steadystate_sensitivity_mode_ {SteadyStateSensitivityMode::newtonOnly};

    /** flag indicating whether reinitialization of states depending on
     *  fixed parameters is activated
     */
    bool reinitialize_fixed_parameter_initial_states_ {false};

    /** Indicates whether the result of every call to Model::f* should be
     * checked for finiteness */
    bool always_check_finite {false};
};

bool operator==(const Model &a, const Model &b);

} // namespace amici

#endif // AMICI_MODEL_H
