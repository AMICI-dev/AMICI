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
void serialize(Archive &ar, amici::Model &m, unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief Exchange format to store and transfer the state of the
 * model at a specific timepoint.
 *
 * This is designed to only encompass the minimal
 * number of attributes that need to be transferred.
 */
struct ModelState {
    /**
     * Flag indicating whether a certain Heaviside function should be active or
     * not (dimension: `ne`)
     */
    std::vector<realtype> h;

    /** Total abundances for conservation laws
     (dimension: `nx_rdata - nx_solver`) */
    std::vector<realtype> total_cl;

    /** Sensitivities of total abundances for conservation laws
     (dimension: `(nx_rdata-nx_solver) x np`, row-major) */
    std::vector<realtype> stotal_cl;

    /** Unscaled parameters (dimension: `np`) */
    std::vector<realtype> unscaledParameters;

    /** Constants (dimension: `nk`) */
    std::vector<realtype> fixedParameters;

    /**
     * Indexes of parameters wrt to which sensitivities are computed
     * (dimension: nplist)
     */
    std::vector<int> plist;
};



/**
 * @brief The Model class represents an AMICI ODE/DAE model.
 *
 * The model can compute various model related quantities based on symbolically
 * generated code.
 */
class Model : public AbstractModel {
  public:
    /** Default constructor */
    Model() = default;

    /**
     * @brief Constructor with model dimensions.
     * @param nx_rdata Number of state variables
     * @param nxtrue_rdata Number of state variables of the non-augmented model
     * @param nx_solver Number of state variables with conservation laws applied
     * @param nxtrue_solver Number of state variables of the non-augmented model
     * with conservation laws applied
     * @param nx_solver_reinit Number of state variables with conservation laws
     * subject to reinitialization
     * @param ny Number of observables
     * @param nytrue Number of observables of the non-augmented model
     * @param nz Number of event observables
     * @param nztrue Number of event observables of the non-augmented model
     * @param ne Number of events
     * @param nJ Number of objective functions
     * @param nw Number of repeating elements
     * @param ndwdx Number of nonzero elements in the `x` derivative of the
     * repeating elements
     * @param ndwdp Number of nonzero elements in the `p` derivative of the
     * repeating elements
     * @param ndwdw Number of nonzero elements in the `w` derivative of the
     * repeating elements
     * @param ndxdotdw Number of nonzero elements in the \f$w\f$ derivative of
     * \f$xdot\f$
     * @param ndJydy Number of nonzero elements in the \f$y\f$ derivative of
     * \f$dJy\f$ (dimension `nytrue`)
     * @param nnz Number of nonzero elements in Jacobian
     * @param ubw Upper matrix bandwidth in the Jacobian
     * @param lbw Lower matrix bandwidth in the Jacobian
     * @param o2mode Second order sensitivity mode
     * @param p Parameters
     * @param k Constants
     * @param plist Indexes wrt to which sensitivities are to be computed
     * @param idlist Indexes indicating algebraic components (DAE only)
     * @param z2event Mapping of event outputs to events
     * @param pythonGenerated Flag indicating matlab or python wrapping
     * @param ndxdotdp_explicit Number of nonzero elements in `dxdotdp_explicit`
     * @param ndxdotdx_explicit Number of nonzero elements in `dxdotdx_explicit`
     * @param w_recursion_depth Recursion depth of fw
     */
    Model(int nx_rdata, int nxtrue_rdata, int nx_solver, int nxtrue_solver,
          int nx_solver_reinit, int ny, int nytrue, int nz, int nztrue, int ne,
          int nJ, int nw, int ndwdx, int ndwdp, int ndwdw, int ndxdotdw,
          std::vector<int> ndJydy, int nnz, int ubw, int lbw,
          amici::SecondOrderMode o2mode,
          const std::vector<amici::realtype> &p, std::vector<amici::realtype> k,
          const std::vector<int> &plist, std::vector<amici::realtype> idlist,
          std::vector<int> z2event, bool pythonGenerated = false,
          int ndxdotdp_explicit = 0, int ndxdotdx_explicit = 0,
          int w_recursion_depth = 0);

    /** Destructor. */
    ~Model() override = default;

    /**
     * @brief Copy assignment is disabled until const members are removed.
     * @param other Object to copy from
     * @return
     */
    Model &operator=(Model const &other) = delete;

    /**
     * @brief Clone this instance.
     * @return The clone
     */
    virtual Model *clone() const = 0;

    /**
     * @brief Serialize Model (see `boost::serialization::serialize`).
     * @param ar Archive to serialize to
     * @param m Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, Model &m,
                                                unsigned int version);

    /**
     * @brief Check equality of data members.
     * @param a First model instance
     * @param b Second model instance
     * @return Equality
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
    using AbstractModel::fdJydy_colptrs;
    using AbstractModel::fdJydy_rowvals;
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
    using AbstractModel::fdwdw;
    using AbstractModel::fdwdw_colptrs;
    using AbstractModel::fdwdw_rowvals;
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
     * @brief Initialize model properties.
     * @param x Reference to state variables
     * @param dx Reference to time derivative of states (DAE only)
     * @param sx Reference to state variable sensitivities
     * @param sdx Reference to time derivative of state sensitivities (DAE only)
     * @param computeSensitivities Flag indicating whether sensitivities are to
     * be computed
     */
    void initialize(AmiVector &x, AmiVector &dx, AmiVectorArray &sx,
                    AmiVectorArray &sdx, bool computeSensitivities);

    /**
     * @brief Initialize model properties.
     * @param xB Adjoint state variables
     * @param dxB Time derivative of adjoint states (DAE only)
     * @param xQB Adjoint quadratures
     * @param posteq Flag indicating whether postequilibration was performed
     */
    void initializeB(AmiVector &xB, AmiVector &dxB, AmiVector &xQB,
                     bool posteq) const;

    /**
     * @brief Initialize initial states.
     * @param x State vector to be initialized
     */
    void initializeStates(AmiVector &x);

    /**
     * @brief Initialize initial state sensitivities.
     * @param sx Reference to state variable sensitivities
     * @param x Reference to state variables
     */
    void initializeStateSensitivities(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Initialize the Heaviside variables `h` at the initial time `t0`.
     *
     * Heaviside variables activate/deactivate on event occurrences.
     *
     * @param x Reference to state variables
     * @param dx Reference to time derivative of states (DAE only)
     */
    void initHeaviside(const AmiVector &x, const AmiVector &dx);

    /**
     * @brief Get number of parameters wrt to which sensitivities are computed.
     * @return Length of sensitivity index vector
     */
    int nplist() const;

    /**
     * @brief Get total number of model parameters.
     * @return Length of parameter vector
     */
    int np() const;

    /**
     * @brief Get number of constants
     * @return Length of constant vector
     */
    int nk() const;

    /**
     * @brief Get number of conservation laws.
     * @return Number of conservation laws (i.e., difference between `nx_rdata`
     * and `nx_solver`).
     */
    int ncl() const;

    /**
     * @brief Get number of solver states subject to reinitialization.
     * @return Model member `nx_solver_reinit`
     */
    int nx_reinit() const;

    /**
     * @brief Get fixed parameters.
     * @return Pointer to constants array
     */
    const double *k() const;

    /**
     * @brief Get maximum number of events that may occur for each type.
     * @return Maximum number of events that may occur for each type
     */
    int nMaxEvent() const;

    /**
     * @brief Set maximum number of events that may occur for each type.
     * @param nmaxevent Maximum number of events that may occur for each type
     */
    void setNMaxEvent(int nmaxevent);

    /**
     * @brief Get number of timepoints.
     * @return Number of timepoints
     */
    int nt() const;

    /**
     * @brief Get parameter scale for each parameter.
     * @return Vector of parameter scales
     */
    std::vector<ParameterScaling> const &getParameterScale() const;

    /**
     * @brief Set parameter scale for each parameter.
     *
     * NOTE: Resets initial state sensitivities.
     *
     * @param pscale Scalar parameter scale to be set for all parameters
     */
    void setParameterScale(ParameterScaling pscale);

    /**
     * @brief Set parameter scale for each parameter.
     *
     * NOTE: Resets initial state sensitivities.
     *
     * @param pscaleVec Vector of parameter scales
     */
    void setParameterScale(const std::vector<ParameterScaling> &pscaleVec);

    /**
     * @brief Get parameters with transformation according to parameter scale
     * applied.
     * @return Unscaled parameters
     */
    std::vector<realtype> const &getUnscaledParameters() const;

    /**
     * @brief Get parameter vector.
     * @return The user-set parameters (see also `Model::getUnscaledParameters`)
     */
    std::vector<realtype> const &getParameters() const;

    /**
     * @brief Get value of first model parameter with the specified ID.
     * @param par_id Parameter ID
     * @return Parameter value
     */
    realtype getParameterById(std::string const &par_id) const;

    /**
     * @brief Get value of first model parameter with the specified name.
     * @param par_name Parameter name
     * @return Parameter value
     */
    realtype getParameterByName(std::string const &par_name) const;

    /**
     * @brief Set the parameter vector.
     * @param p Vector of parameters
     */
    void setParameters(std::vector<realtype> const &p);

    /**
     * @brief Set model parameters according to the parameter IDs and mapped
     * values.
     * @param p Map of parameters IDs and values
     * @param ignoreErrors Ignore errors such as parameter IDs in p which are
     * not model parameters
     */
    void setParameterById(std::map<std::string, realtype> const &p,
                          bool ignoreErrors = false);

    /**
     * @brief Set value of first model parameter with the specified ID.
     * @param par_id Parameter ID
     * @param value Parameter value
     */
    void setParameterById(std::string const &par_id, realtype value);

    /**
     * @brief Set all values of model parameters with IDs matching the specified
     * regular expression.
     * @param par_id_regex Parameter ID regex
     * @param value Parameter value
     * @return Number of parameter IDs that matched the regex
     */
    int setParametersByIdRegex(std::string const &par_id_regex, realtype value);

    /**
     * @brief Set value of first model parameter with the specified name.
     * @param par_name Parameter name
     * @param value Parameter value
     */
    void setParameterByName(std::string const &par_name, realtype value);

    /**
     * @brief Set model parameters according to the parameter name and mapped
     * values.
     * @param p Map of parameters names and values
     * @param ignoreErrors Ignore errors such as parameter names in p which are
     * not model parameters
     */
    void setParameterByName(std::map<std::string, realtype> const &p,
                            bool ignoreErrors = false);

    /**
     * @brief Set all values of all model parameters with names matching the
     * specified regex.
     * @param par_name_regex Parameter name regex
     * @param value Parameter value
     * @return Number of fixed parameter names that matched the regex
     */
    int setParametersByNameRegex(std::string const &par_name_regex,
                                 realtype value);

    /**
     * @brief Get values of fixed parameters.
     * @return Vector of fixed parameters with same ordering as in
     * Model::getFixedParameterIds
     */
    std::vector<realtype> const &getFixedParameters() const;

    /**
     * @brief Get value of fixed parameter with the specified ID.
     * @param par_id Parameter ID
     * @return Parameter value
     */
    realtype getFixedParameterById(std::string const &par_id) const;

    /**
     * @brief Get value of fixed parameter with the specified name.
     *
     * If multiple parameters have the same name, the first parameter with
     * matching name is returned.
     *
     * @param par_name Parameter name
     * @return Parameter value
     */
    realtype getFixedParameterByName(std::string const &par_name) const;

    /**
     * @brief Set values for constants.
     * @param k Vector of fixed parameters
     */
    void setFixedParameters(std::vector<realtype> const &k);

    /**
     * @brief Set value of first fixed parameter with the specified ID.
     * @param par_id Fixed parameter id
     * @param value Fixed parameter value
     */
    void setFixedParameterById(std::string const &par_id, realtype value);

    /**
     * @brief Set values of all fixed parameters with the ID matching the
     * specified regex.
     * @param par_id_regex Fixed parameter name regex
     * @param value Fixed parameter value
     * @return Number of fixed parameter IDs that matched the regex
     */
    int setFixedParametersByIdRegex(std::string const &par_id_regex,
                                    realtype value);

    /**
     * @brief Set value of first fixed parameter with the specified name.
     * @param par_name Fixed parameter ID
     * @param value Fixed parameter value
     */
    void setFixedParameterByName(std::string const &par_name, realtype value);

    /**
     * @brief Set value of all fixed parameters with name matching the specified
     * regex.
     * @param par_name_regex Fixed parameter name regex
     * @param value Fixed parameter value
     * @return Number of fixed parameter names that matched the regex
     */
    int setFixedParametersByNameRegex(std::string const &par_name_regex,
                                      realtype value);

    /**
     * @brief Get the model name.
     * @return Model name
     */
    virtual std::string getName() const;

    /**
     * @brief Report whether the model has parameter names set.
     *
     * @return Boolean indicating whether parameter names were set. Also returns
     * `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasParameterNames() const;

    /**
     * @brief Get names of the model parameters.
     * @return The parameter names
     */
    virtual std::vector<std::string> getParameterNames() const;

    /**
     * @brief Report whether the model has state names set.
     *
     * @return Boolean indicating whether state names were set. Also returns
     * `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasStateNames() const;

    /**
     * @brief Get names of the model states.
     * @return State names
     */
    virtual std::vector<std::string> getStateNames() const;

    /**
     * @brief Report whether the model has fixed parameter names set.
     * @return Boolean indicating whether fixed parameter names were set. Also
     * returns `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasFixedParameterNames() const;

    /**
     * @brief Get names of the fixed model parameters.
     * @return Fixed parameter names
     */
    virtual std::vector<std::string> getFixedParameterNames() const;

    /**
     * @brief Report whether the model has observable names set.
     * @return Boolean indicating whether observable names were set. Also
     * returns `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasObservableNames() const;

    /**
     * @brief Get names of the observables.
     * @return Observable names
     */
    virtual std::vector<std::string> getObservableNames() const;

    /**
     * @brief Report whether the model has parameter IDs set.
     * @return Boolean indicating whether parameter IDs were set. Also returns
     * `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasParameterIds() const;

    /**
     * @brief Get IDs of the model parameters.
     * @return Parameter IDs
     */
    virtual std::vector<std::string> getParameterIds() const;

    /**
     * @brief Report whether the model has state IDs set.
     * @return Boolean indicating whether state IDs were set. Also returns
     * `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasStateIds() const;

    /**
     * @brief Get IDs of the model states.
     * @return Sate IDs
     */
    virtual std::vector<std::string> getStateIds() const;

    /**
     * @brief Report whether the model has fixed parameter IDs set.
     * @return Boolean indicating whether fixed parameter IDs were set. Also
     * returns `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasFixedParameterIds() const;

    /**
     * @brief Get IDs of the fixed model parameters.
     * @return Fixed parameter IDs
     */
    virtual std::vector<std::string> getFixedParameterIds() const;

    /**
     * @brief Report whether the model has observable IDs set.
     * @return Boolean indicating whether observable ids were set. Also returns
     * `true` if the number of corresponding variables is just zero.
     */
    virtual bool hasObservableIds() const;

    /**
     * @brief Get IDs of the observables.
     * @return Observable IDs
     */
    virtual std::vector<std::string> getObservableIds() const;
    
    /**
     * @brief Checks whether the defined noise model is gaussian, i.e., the nllh is quadratic
     * @return boolean flag
     */
    virtual bool hasQuadraticLLH() const;

    /**
     * @brief Get the timepoint vector.
     * @return Timepoint vector
     */
    std::vector<realtype> const &getTimepoints() const;

    /**
     * @brief Get simulation timepoint for time index `it`.
     * @param it Time index
     * @return Timepoint
     */
    realtype getTimepoint(int it) const;

    /**
     * @brief Set the timepoint vector.
     * @param ts New timepoint vector
     */
    void setTimepoints(std::vector<realtype> const &ts);

    /**
     * @brief Get simulation start time.
     * @return Simulation start time
     */
    double t0() const;

    /**
     * @brief Set simulation start time.
     * @param t0 Simulation start time
     */
    void setT0(double t0);

    /**
     * @brief Get flags indicating whether states should be treated as
     * non-negative.
     * @return Vector of flags
     */
    std::vector<bool> const &getStateIsNonNegative() const;

    /**
     * @brief Set flags indicating whether states should be treated as
     * non-negative.
     * @param stateIsNonNegative Vector of flags
     */
    void setStateIsNonNegative(std::vector<bool> const &stateIsNonNegative);

    /**
     * @brief Set flags indicating that all states should be treated as
     * non-negative.
     */
    void setAllStatesNonNegative();

    /**
     * @brief Get the current model state.
     * @return Current model state
     */
    ModelState const &getModelState() const {
        return state_;
    };

    /**
     * @brief Set the current model state.
     * @param state Model state
     */
    void setModelState(ModelState const &state) {
        if (static_cast<int>(state.unscaledParameters.size()) != np())
            throw AmiException("Mismatch in parameter size");
        if (static_cast<int>(state.fixedParameters.size()) != nk())
            throw AmiException("Mismatch in fixed parameter size");
        if (static_cast<int>(state.h.size()) != ne)
            throw AmiException("Mismatch in Heaviside size");
        if (static_cast<int>(state.total_cl.size()) != ncl())
            throw AmiException("Mismatch in conservation law size");
        if (static_cast<int>(state.stotal_cl.size()) != ncl() * np() )
            throw AmiException("Mismatch in conservation law sensitivity size");
        state_ = state;
    };

    /**
     * @brief Get the list of parameters for which sensitivities are computed.
     * @return List of parameter indices
     */
    std::vector<int> const &getParameterList() const;

    /**
     * @brief Get entry in parameter list by index.
     * @param pos Index in sensitivity parameter list
     * @return Index in parameter list
     */
    int plist(int pos) const;

    /**
     * @brief Set the list of parameters for which sensitivities are to be
     * computed.
     *
     * NOTE: Resets initial state sensitivities.
     *
     * @param plist List of parameter indices
     */
    void setParameterList(std::vector<int> const &plist);

    /**
     * @brief Get the initial states.
     * @return Initial state vector
     */
    std::vector<realtype> getInitialStates();

    /**
     * @brief Set the initial states.
     * @param x0 Initial state vector
     */
    void setInitialStates(std::vector<realtype> const &x0);

    /**
     * @brief Return whether custom initial states have been set.
     * @return `true` if has custom initial states, otherwise `false`
     */
    bool hasCustomInitialStates() const;

    /**
     * @brief Get the initial states sensitivities.
     * @return vector of initial state sensitivities
     */
    std::vector<realtype> getInitialStateSensitivities();

    /**
     * @brief Set the initial state sensitivities.
     * @param sx0 vector of initial state sensitivities with chainrule applied.
     * This could be a slice of ReturnData::sx or ReturnData::sx0
     */
    void setInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief Return whether custom initial state sensitivities have been set.
     * @return `true` if has custom initial state sensitivities, otherwise
     * `false`.
     */
    bool hasCustomInitialStateSensitivities() const;

    /**
     * @brief Set the initial state sensitivities.
     * @param sx0 Vector of initial state sensitivities without chainrule
     * applied. This could be the readin from a `model.sx0data` saved to HDF5.
     */
    void setUnscaledInitialStateSensitivities(std::vector<realtype> const &sx0);

    /**
     * @brief Set the mode how sensitivities are computed in the steadystate
     * simulation.
     * @param mode Steadystate sensitivity mode
     */
    void setSteadyStateSensitivityMode(SteadyStateSensitivityMode mode);

    /**
     * @brief Gets the mode how sensitivities are computed in the steadystate
     * simulation.
     * @return Mode
     */
    SteadyStateSensitivityMode getSteadyStateSensitivityMode() const;

    /**
     * @brief Set whether initial states depending on fixed parameters are to be
     * reinitialized after preequilibration and presimulation.
     * @param flag Fixed parameters reinitialized?
     */
    void setReinitializeFixedParameterInitialStates(bool flag);

    /**
     * @brief Get whether initial states depending on fixedParameters are to be
     * reinitialized after preequilibration and presimulation.
     * @return flag `true` / `false`
     */
    bool getReinitializeFixedParameterInitialStates() const;

    /**
     * @brief Require computation of sensitivities for all parameters p [0..np[
     * in natural order.
     *
     * NOTE: Resets initial state sensitivities.
     */
    void requireSensitivitiesForAllParameters();

    /**
     * @brief Get time-resolved `w`.
     * @param w Buffer (dimension: `nw`)
     * @param t Current timepoint
     * @param x Current state
     */
    void getExpression(gsl::span<realtype> w, const realtype t, const AmiVector &x);

    /**
     * @brief Get time-resolved observables.
     * @param y Buffer (dimension: `ny`)
     * @param t Current timepoint
     * @param x Current state
     */
    void getObservable(gsl::span<realtype> y, const realtype t,
                       const AmiVector &x);

    /**
     * @brief Get sensitivity of time-resolved observables.
     *
     * Total derivative \f$sy = dydx * sx + dydp\f$
     * (only for forward sensitivities).
     * @param sy buffer (dimension: `ny` x `nplist`, row-major)
     * @param t Timpoint
     * @param x State variables
     * @param sx State sensitivities
     */
    void getObservableSensitivity(gsl::span<realtype> sy, const realtype t,
                                  const AmiVector &x, const AmiVectorArray &sx);

    /**
     * @brief Get time-resolved observable standard deviations
     * @param sigmay Buffer (dimension: `ny`)
     * @param it Timepoint index
     * @param edata Pointer to experimental data instance (optional, pass
     * `nullptr` to ignore)
     */
    void getObservableSigma(gsl::span<realtype> sigmay, const int it,
                            const ExpData *edata);

    /**
     * @brief Sensitivity of time-resolved observable standard deviation.
     *
     * Total derivative (can be used with both adjoint and forward sensitivity).
     *
     * @param ssigmay Buffer (dimension: `ny` x `nplist`, row-major)
     * @param it Timepoint index
     * @param edata Pointer to experimental data instance (optional, pass
     * `nullptr` to ignore)
     */
    void getObservableSigmaSensitivity(gsl::span<realtype> ssigmay,
                                       const int it, const ExpData *edata);

    /**
     * @brief Add time-resolved measurement negative log-likelihood \f$Jy\f$.
     * @param Jy Buffer (dimension: 1)
     * @param it Timepoint index
     * @param x State variables
     * @param edata Experimental data
     */
    void addObservableObjective(realtype &Jy, const int it, const AmiVector &x,
                                const ExpData &edata);

    /**
     * @brief Add sensitivity of time-resolved measurement negative log-likelihood
     * \f$Jy\f$.
     *
     * @param sllh First-order buffer (dimension: `nplist`)
     * @param s2llh Second-order buffer (dimension: `nJ - 1` x `nplist`, row-major)
     * @param it Timepoint index
     * @param x State variables
     * @param sx State sensitivities
     * @param edata Experimental data
     */
    void addObservableObjectiveSensitivity(std::vector<realtype> &sllh,
                                           std::vector<realtype> &s2llh,
                                           const int it, const AmiVector &x,
                                           const AmiVectorArray &sx,
                                           const ExpData &edata);

    /**
     * @brief Add sensitivity of time-resolved measurement negative
     * log-likelihood \f$Jy\f$.
     *
     * Partial derivative (to be used with adjoint sensitivities).
     *
     * @param sllh First order output buffer (dimension: `nplist`)
     * @param s2llh Second order output buffer
     * (dimension: `nJ - 1` x `nplist`, row-major)
     * @param it Timepoint index
     * @param x State variables
     * @param edata Experimental data
     */
    void addPartialObservableObjectiveSensitivity(std::vector<realtype> &sllh,
                                                  std::vector<realtype> &s2llh,
                                                  const int it,
                                                  const AmiVector &x,
                                                  const ExpData &edata);

    /**
     * @brief Get state sensitivity of the negative loglikelihood \f$Jy\f$,
     * partial derivative (to be used with adjoint sensitivities).
     *
     * @param dJydx Output buffer (dimension: `nJ` x `nx_solver`, row-major)
     * @param it Timepoint index
     * @param x State variables
     * @param edata Experimental data instance
     */
    void getAdjointStateObservableUpdate(gsl::span<realtype> dJydx,
                                         const int it, const AmiVector &x,
                                         const ExpData &edata);

    /**
     * @brief Get event-resolved observables.
     * @param z Output buffer (dimension: `nz`)
     * @param ie Event index
     * @param t Timepoint
     * @param x State variables
     */
    void getEvent(gsl::span<realtype> z, const int ie, const realtype t,
                  const AmiVector &x);
    /**
     * @brief Get sensitivities of event-resolved observables.
     *
     * Total derivative (only forward sensitivities).
     *
     * @param sz Output buffer (dimension: `nz x nplist`, row-major)
     * @param ie Event index
     * @param t Timepoint
     * @param x State variables
     * @param sx State sensitivities
     */
    void getEventSensitivity(gsl::span<realtype> sz, const int ie,
                             const realtype t, const AmiVector &x,
                             const AmiVectorArray &sx);

    /**
     * @brief Get sensitivity of `z` at final timepoint.
     *
     * Ignores sensitivity of timepoint. Total derivative.
     *
     * @param sz Output buffer (dimension: `nz x nplist`, row-major)
     * @param ie Event index
     */
    void getUnobservedEventSensitivity(gsl::span<realtype> sz, const int ie);

    /**
     * @brief Get regularization for event-resolved observables.
     * @param rz Output buffer (dimension: `nz`)
     * @param ie Event index
     * @param t Timepoint
     * @param x State variables
     */
    void getEventRegularization(gsl::span<realtype> rz, const int ie,
                                const realtype t, const AmiVector &x);

    /**
     * @brief Get sensitivities of regularization for event-resolved
     * observables.
     *
     * Total derivative. Only forward sensitivities.
     *
     * @param srz Output buffer (dimension: `nz x nplist`, row-major)
     * @param ie Event index
     * @param t Timepoint
     * @param x State variables
     * @param sx State sensitivities
     */
    void getEventRegularizationSensitivity(gsl::span<realtype> srz,
                                           const int ie, const realtype t,
                                           const AmiVector &x,
                                           const AmiVectorArray &sx);
    /**
     * @brief Get event-resolved observable standard deviations.
     * @param sigmaz Output buffer (dimension: `nz`)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param edata Pointer to experimental data (optional, pass
     * `nullptr` to ignore)
     */
    void getEventSigma(gsl::span<realtype> sigmaz, const int ie,
                       const int nroots, const realtype t,
                       const ExpData *edata);

    /**
     * @brief Get sensitivities of event-resolved observable standard
     * deviations.
     *
     * Total derivative (only forward sensitivities).
     *
     * @param ssigmaz Output buffer (dimension: `nz x nplist`, row-major)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param edata Pointer to experimental data (optional, pass
     * `nullptr` to ignore)
     */
    void getEventSigmaSensitivity(gsl::span<realtype> ssigmaz, const int ie,
                                  const int nroots, const realtype t,
                                  const ExpData *edata);

    /**
     * @brief Add event-resolved observable negative log-likelihood.
     * @param Jz Output buffer (dimension: 1)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void addEventObjective(realtype &Jz, const int ie, const int nroots,
                           const realtype t, const AmiVector &x,
                           const ExpData &edata);

    /**
     * @brief Add event-resolved observable negative log-likelihood.
     * @param Jrz Output buffer (dimension: 1)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void addEventObjectiveRegularization(realtype &Jrz, const int ie,
                                         const int nroots, const realtype t,
                                         const AmiVector &x,
                                         const ExpData &edata);

    /**
     * @brief Add sensitivity of time-resolved measurement negative
     * log-likelihood \f$Jy\f$.
     *
     * Total derivative (to be used with forward sensitivities).
     *
     * @param sllh First order buffer (dimension: `nplist`)
     * @param s2llh Second order buffer
     * (dimension: `(nJ-1) x nplist`, row-major)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param x State variables
     * @param sx State sensitivities
     * @param edata Experimental data
     */
    void addEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                      std::vector<realtype> &s2llh,
                                      const int ie, const int nroots,
                                      const realtype t, const AmiVector &x,
                                      const AmiVectorArray &sx,
                                      const ExpData &edata);

    /**
     * @brief Add sensitivity of time-resolved measurement negative
     * log-likelihood \f$Jy\f$.
     *
     * Partial derivative (to be used with adjoint sensitivities).
     *
     * @param sllh First order buffer (dimension: `nplist`)
     * @param s2llh Second order buffer
     * (dimension: `(nJ-1) x nplist`, row-major)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void addPartialEventObjectiveSensitivity(std::vector<realtype> &sllh,
                                             std::vector<realtype> &s2llh,
                                             const int ie, const int nroots,
                                             const realtype t,
                                             const AmiVector &x,
                                             const ExpData &edata);

    /**
     * @brief State sensitivity of the negative loglikelihood \f$Jz\f$.
     *
     * Partial derivative (to be used with adjoint sensitivities).
     *
     * @param dJzdx Output buffer (dimension: `nJ x nx_solver`, row-major)
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void getAdjointStateEventUpdate(gsl::span<realtype> dJzdx, const int ie,
                                    const int nroots, const realtype t,
                                    const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of event timepoint, total derivative.
     *
     * Only forward sensitivities.
     *
     * @param stau Timepoint sensitivity (dimension: `nplist`)
     * @param t Timepoint
     * @param ie Event index
     * @param x State variables
     * @param sx State sensitivities
     */
    void getEventTimeSensitivity(std::vector<realtype> &stau, const realtype t,
                                 const int ie, const AmiVector &x,
                                 const AmiVectorArray &sx);

    /**
     * @brief Update state variables after event.
     * @param x Current state (will be overwritten)
     * @param ie Event index
     * @param t Current timepoint
     * @param xdot Current residual function values
     * @param xdot_old Value of residual function before event
     */
    void addStateEventUpdate(AmiVector &x, const int ie, const realtype t,
                             const AmiVector &xdot, const AmiVector &xdot_old);

    /**
     * @brief Update state sensitivity after event.
     * @param sx Current state sensitivity (will be overwritten)
     * @param ie Event index
     * @param t Current timepoint
     * @param x_old Current state
     * @param xdot Current residual function values
     * @param xdot_old Value of residual function before event
     * @param stau Timepoint sensitivity, to be computed with
     * `Model::getEventTimeSensitivity`
     */
    void addStateSensitivityEventUpdate(AmiVectorArray &sx, const int ie,
                                        const realtype t,
                                        const AmiVector &x_old,
                                        const AmiVector &xdot,
                                        const AmiVector &xdot_old,
                                        const std::vector<realtype> &stau);

    /**
     * @brief Update adjoint state after event.
     * @param xB Current adjoint state (will be overwritten)
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     * @param xdot Current residual function values
     * @param xdot_old Value of residual function before event
     */
    void addAdjointStateEventUpdate(AmiVector &xB, const int ie,
                                    const realtype t, const AmiVector &x,
                                    const AmiVector &xdot,
                                    const AmiVector &xdot_old);

    /**
     * @brief Update adjoint quadratures after event.
     * @param xQB Current quadrature state (will be overwritten)
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     * @param xB Current adjoint state
     * @param xdot Current residual function values
     * @param xdot_old Value of residual function before event
     */
    void addAdjointQuadratureEventUpdate(AmiVector xQB, const int ie,
                                         const realtype t, const AmiVector &x,
                                         const AmiVector &xB,
                                         const AmiVector &xdot,
                                         const AmiVector &xdot_old);

    /**
     * @brief Update the Heaviside variables `h` on event occurrences.
     *
     * @param rootsfound Provides the direction of the zero-crossing, so adding
     * it will give the right update to the Heaviside variables (zero if no root
     * was found)
     */
    void updateHeaviside(const std::vector<int> &rootsfound);

    /**
     * @brief Updates the Heaviside variables `h` on event occurrences in the
     * backward problem.
     * @param rootsfound Provides the direction of the zero-crossing, so adding
     * it will give the right update to the Heaviside variables (zero if no root
     * was found)
     */
    void updateHeavisideB(const int *rootsfound);

    /**
     * @brief Check if the given array has only finite elements.
     *
     * If not, try to give hints by which other fields this could be caused.
     *
     * @param array Array to check
     * @param fun Name of the function that generated the values (for more
     * informative messages).
     * @return `amici::AMICI_RECOVERABLE_ERROR` if a NaN/Inf value was found,
     * `amici::AMICI_SUCCESS` otherwise
     */
    int checkFinite(gsl::span<const realtype> array, const char *fun) const;

    /**
     * @brief Set whether the result of every call to `Model::f*` should be
     * checked for finiteness.
     * @param alwaysCheck
     */
    void setAlwaysCheckFinite(bool alwaysCheck);

    /**
     * @brief Get setting of whether the result of every call to `Model::f*`
     * should be checked for finiteness.
     * @return that
     */
    bool getAlwaysCheckFinite() const;

    /**
     * @brief Compute/get initial states.
     * @param x Output buffer.
     */
    void fx0(AmiVector &x);

    /**
     * @brief Set only those initial states that are specified via
     * fixed parameters.
     * @param x Output buffer.
     */
    void fx0_fixedParameters(AmiVector &x);

    /**
     * @brief Compute/get initial value for initial state sensitivities.
     * @param sx Output buffer for state sensitivities
     * @param x State variables
     */
    void fsx0(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Get only those initial states sensitivities that are affected
     * from `amici::Model::fx0_fixedParameters`.
     * @param sx Output buffer for state sensitivities
     * @param x State variables
     */
    void fsx0_fixedParameters(AmiVectorArray &sx, const AmiVector &x);

    /**
     * @brief Compute sensitivity of derivative initial states sensitivities
     * `sdx0`.
     *
     * Only necessary for DAEs.
     */
    virtual void fsdx0();

    /**
     * @brief Expand conservation law for states.
     * @param x_rdata Output buffer for state variables with conservation laws
     * expanded (stored in `amici::ReturnData`).
     * @param x_solver State variables with conservation laws applied
     * (solver returns this)
     */
    void fx_rdata(AmiVector &x_rdata, const AmiVector &x_solver);

    /**
     * @brief Expand conservation law for state sensitivities.
     * @param sx_rdata Output buffer for state variables sensitivities with
     * conservation laws expanded (stored in `amici::ReturnData`).
     * @param sx_solver State variables sensitivities with conservation laws
     * applied (solver returns this)
     */
    void fsx_rdata(AmiVectorArray &sx_rdata, const AmiVectorArray &sx_solver);

    /** Number of states */
    int nx_rdata{0};

    /** Number of states in the unaugmented system */
    int nxtrue_rdata{0};

    /** Number of states with conservation laws applied */
    int nx_solver{0};

    /**
     * Number of states in the unaugmented system with conservation laws
     * applied
     */
    int nxtrue_solver{0};

    /** Number of solver states subject to reinitialization */
    int nx_solver_reinit{0};

    /** Number of observables */
    int ny{0};

    /** Number of observables in the unaugmented system */
    int nytrue{0};

    /** Number of event outputs */
    int nz{0};

    /** Number of event outputs in the unaugmented system */
    int nztrue{0};

    /** Number of events */
    int ne{0};

    /** Number of common expressions */
    int nw{0};

    /** Number of nonzero entries in Jacobian */
    int nnz{0};

    /** Dimension of the augmented objective function for 2nd order ASA */
    int nJ{0};

    /** Upper bandwidth of the Jacobian */
    int ubw{0};

    /** Lower bandwidth of the Jacobian */
    int lbw{0};

    /** Flag indicating Matlab- or Python-based model generation */
    bool pythonGenerated;
    
    /**
     * @brief getter for dxdotdp (matlab generated)
     * @return dxdotdp
     */
    const AmiVectorArray &get_dxdotdp() const;
    
    /**
     * @brief getter for dxdotdp (python generated)
     * @return dxdotdp
     */
    const SUNMatrixWrapper &get_dxdotdp_full() const;

    /**
     * Flag indicating whether for
     * `amici::Solver::sensi_` == `amici::SensitivityOrder::second`
     * directional or full second order derivative will be computed
     */
    SecondOrderMode o2mode{SecondOrderMode::none};

    /** Flag array for DAE equations */
    std::vector<realtype> idlist;

    /** AMICI application context */
    AmiciApplication *app = &defaultContext;

  protected:
    /**
     * @brief Write part of a slice to a buffer according to indices specified
     * in z2event.
     * @param slice Input data slice
     * @param buffer Output data slice
     * @param ie Event index
     */
    void writeSliceEvent(gsl::span<const realtype> slice,
                         gsl::span<realtype> buffer, const int ie);

    /**
     * @brief Write part of a sensitivity slice to a buffer according to
     * indices specified in z2event.
     * @param slice source data slice
     * @param buffer output data slice
     * @param ie event index
     */
    void writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                    gsl::span<realtype> buffer, const int ie);

    /**
     * @brief Separate first and second order objective sensitivity information
     * and write them into the respective buffers.
     * @param dLLhdp Data with mangled first- and second-order information
     * @param sllh First order buffer
     * @param s2llh Second order buffer
     */
    void writeLLHSensitivitySlice(const std::vector<realtype> &dLLhdp,
                                  std::vector<realtype> &sllh,
                                  std::vector<realtype> &s2llh);

    /**
     * @brief Verify that the provided buffers have the expected size.
     * @param sllh first order buffer
     * @param s2llh second order buffer
     */
    void checkLLHBufferSize(const std::vector<realtype> &sllh,
                            const std::vector<realtype> &s2llh) const;

    /**
     * @brief Set the nplist-dependent vectors to their proper sizes.
     */
    void initializeVectors();

    /**
     * @brief Compute observables / measurements.
     * @param t Current timepoint
     * @param x Current state
     */
    void fy(realtype t, const AmiVector &x);

    /**
     * @brief Compute partial derivative of observables \f$y\f$ w.r.t. model
     * parameters `p`.
     * @param t Current timepoint
     * @param x Current state
     */
    void fdydp(realtype t, const AmiVector &x);

    /**
     * @brief Compute partial derivative of observables \f$y\f$ w.r.t. state
     * variables `x`.
     * @param t Current timepoint
     * @param x Current state
     */
    void fdydx(realtype t, const AmiVector &x);

    /**
     * @brief Compute standard deviation of measurements.
     * @param it Timepoint index
     * @param edata Experimental data
     */
    void fsigmay(int it, const ExpData *edata);

    /**
     * @brief Compute partial derivative of standard deviation of measurements
     * w.r.t. model parameters.
     * @param it Timepoint index
     * @param edata pointer to `amici::ExpData` data instance holding sigma values
     */
    void fdsigmaydp(int it, const ExpData *edata);

    /**
     * @brief Compute negative log-likelihood of measurements \f$y\f$.
     *
     * @param Jy Variable to which llh will be added
     * @param it Timepoint index
     * @param y Simulated observable
     * @param edata Pointer to experimental data instance
     */
    void fJy(realtype &Jy, int it, const AmiVector &y, const ExpData &edata);

    /**
     * @brief Compute partial derivative of time-resolved measurement negative
     * log-likelihood \f$Jy\f$.
     * @param it timepoint index
     * @param x state variables
     * @param edata Pointer to experimental data
     */
    void fdJydy(int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * Jy w.r.t. standard deviation sigma.
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydsigma(int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute sensitivity of time-resolved measurement negative
     * log-likelihood \f$Jy\f$ w.r.t. parameters for the given timepoint.
     * @param it timepoint index
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJydp(const int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Sensitivity of time-resolved measurement negative log-likelihood
     * \f$Jy\f$ w.r.t. state variables.
     * @param it Timepoint index
     * @param x State variables
     * @param edata Pointer to experimental data instance
     */
    void fdJydx(const int it, const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute event-resolved output.
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     */
    void fz(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute partial derivative of event-resolved output `z` w.r.t.
     * model parameters `p`
     * @param ie event index
     * @param t current timepoint
     * @param x current state
     */
    void fdzdp(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute partial derivative of event-resolved output `z` w.r.t.
     * model states `x`.
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     */
    void fdzdx(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute event root function of events.
     *
     * Equal to `Model::froot` but does not include non-output events.
     *
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     */
    void frz(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute sensitivity of event-resolved root output w.r.t. model
     * parameters `p`.
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     */
    void fdrzdp(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute sensitivity of event-resolved measurements \f$rz\f$ w.r.t.
     * model states `x`.
     * @param ie Event index
     * @param t Current timepoint
     * @param x Current state
     */
    void fdrzdx(int ie, realtype t, const AmiVector &x);

    /**
     * @brief Compute standard deviation of events.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param edata Experimental data
     */
    void fsigmaz(const int ie, const int nroots, const realtype t,
                 const ExpData *edata);

    /**
     * @brief Compute sensitivity of standard deviation of events measurements
     * w.r.t. model parameters `p`.
     * @param ie Event index
     * @param nroots Event occurrence
     * @param t Current timepoint
     * @param edata Pointer to experimental data instance
     */
    void fdsigmazdp(int ie, int nroots, realtype t, const ExpData *edata);

    /**
     * @brief Compute negative log-likelihood of event-resolved measurements
     * `z`.
     * @param Jz Variable to which llh will be added
     * @param nroots Event index
     * @param z Simulated event
     * @param edata Experimental data
     */
    void fJz(realtype &Jz, int nroots, const AmiVector &z,
             const ExpData &edata);

    /**
     * @brief Compute partial derivative of event measurement negative
     * log-likelihood \f$Jz\f$.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void fdJzdz(const int ie, const int nroots, const realtype t,
                const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute sensitivity of event measurement negative log-likelihood
     * \f$Jz\f$ w.r.t. standard deviation sigmaz.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param x State variables
     * @param edata Pointer to experimental data instance
     */
    void fdJzdsigma(const int ie, const int nroots, const realtype t,
                    const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute sensitivity of event-resolved measurement negative
     * log-likelihood Jz w.r.t. parameters.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param x State variables
     * @param edata Pointer to experimental data instance
     */
    void fdJzdp(const int ie, const int nroots, realtype t, const AmiVector &x,
                const ExpData &edata);

    /**
     * @brief Compute sensitivity of event-resolved measurement negative
     * log-likelihood Jz w.r.t. state variables.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void fdJzdx(const int ie, const int nroots, realtype t, const AmiVector &x,
                const ExpData &edata);

    /**
     * @brief Compute regularization of negative log-likelihood with roots of
     * event-resolved measurements rz.
     * @param Jrz Variable to which regularization will be added
     * @param nroots Event index
     * @param rz Regularization variable
     * @param edata Experimental data
     */
    void fJrz(realtype &Jrz, int nroots, const AmiVector &rz,
              const ExpData &edata);

    /**
     * @brief Compute partial derivative of event measurement negative
     * log-likelihood J.
     * @param ie Event index
     * @param nroots Event index
     * @param t Current timepoint
     * @param x State variables
     * @param edata Experimental data
     */
    void fdJrzdz(const int ie, const int nroots, const realtype t,
                 const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute sensitivity of event measurement negative log-likelihood
     * Jz w.r.t. standard deviation sigmaz
     * @param ie event index
     * @param nroots event index
     * @param t current timepoint
     * @param x state variables
     * @param edata pointer to experimental data instance
     */
    void fdJrzdsigma(const int ie, const int nroots, const realtype t,
                     const AmiVector &x, const ExpData &edata);

    /**
     * @brief Compute recurring terms in xdot.
     * @param t Timepoint
     * @param x Array with the states
     */
    void fw(realtype t, const realtype *x);

    /**
     * @brief Compute parameter derivative for recurring terms in xdot.
     * @param t Timepoint
     * @param x Array with the states
     */
    void fdwdp(realtype t, const realtype *x);

    /**
     * @brief Compute state derivative for recurring terms in xdot.
     * @param t Timepoint
     * @param x Array with the states
     */
    void fdwdx(realtype t, const realtype *x);
    
    /**
     * @brief Compute self derivative for recurring terms in xdot.
     * @param t Timepoint
     * @param x Array with the states
     */
    void fdwdw(realtype t, const realtype *x);

    /**
     * @brief Compute fx_rdata.
     *
     * To be implemented by derived class if applicable.
     *
     * @param x_rdata State variables with conservation laws expanded
     * @param x_solver State variables with conservation laws applied
     * @param tcl Total abundances for conservation laws
     */
    virtual void fx_rdata(realtype *x_rdata, const realtype *x_solver,
                          const realtype *tcl);

    /**
     * @brief Compute fsx_solver.
     *
     * To be implemented by derived class if applicable.
     *
     * @param sx_rdata State sensitivity variables with conservation laws
     * expanded
     * @param sx_solver State sensitivity variables with conservation laws
     * applied
     * @param stcl Sensitivities of total abundances for conservation laws
     * @param ip Sensitivity index
     */
    virtual void fsx_rdata(realtype *sx_rdata, const realtype *sx_solver,
                           const realtype *stcl, int ip);

    /**
     * @brief Compute fx_solver.
     *
     * To be implemented by derived class if applicable.
     *
     * @param x_solver State variables with conservation laws applied
     * @param x_rdata State variables with conservation laws expanded
     */
    virtual void fx_solver(realtype *x_solver, const realtype *x_rdata);

    /**
     * @brief Compute fsx_solver.
     *
     * To be implemented by derived class if applicable.
     *
     * @param sx_rdata State sensitivity variables with conservation laws
     * expanded
     * @param sx_solver State sensitivity variables with conservation laws
     * applied
     */
    virtual void fsx_solver(realtype *sx_solver, const realtype *sx_rdata);

    /**
     * @brief Compute ftotal_cl.
     *
     * To be implemented by derived class if applicable.
     *
     * @param total_cl Total abundances of conservation laws
     * @param x_rdata State variables with conservation laws expanded
     */
    virtual void ftotal_cl(realtype *total_cl, const realtype *x_rdata);

    /**
     * @brief Compute fstotal_cl
     *
     * To be implemented by derived class if applicable.
     *
     * @param stotal_cl Sensitivities for the total abundances of conservation
     * laws
     * @param sx_rdata State sensitivity variables with conservation laws
     * expanded
     * @param ip Sensitivity index
     */
    virtual void fstotal_cl(realtype *stotal_cl, const realtype *sx_rdata,
                            int ip);

    /**
     * @brief Compute non-negative state vector.
     *
     * Compute non-negative state vector according to stateIsNonNegative.
     * If anyStateNonNegative is set to `false`, i.e., all entries in
     * stateIsNonNegative are `false`, this function directly returns `x`,
     * otherwise all entries of x are copied in to `amici::Model::x_pos_tmp_`
     * and negative values are replaced by `0` where applicable.
     *
     * @param x State vector possibly containing negative values
     * @return State vector with negative values replaced by `0` according to
     * stateIsNonNegative
     */
    N_Vector computeX_pos(const_N_Vector x);

    /** All variables necessary for function evaluation */
    ModelState state_;

    /** Sparse Jacobian (dimension: `amici::Model::nnz`) */
    mutable SUNMatrixWrapper J_;

    /** Sparse Backwards Jacobian (dimension: `amici::Model::nnz`) */
    mutable SUNMatrixWrapper JB_;

    /** Sparse dxdotdw temporary storage (dimension: `ndxdotdw`) */
    mutable SUNMatrixWrapper dxdotdw_;

    /** Sparse dwdx temporary storage (dimension: `ndwdx`) */
    mutable SUNMatrixWrapper dwdx_;
    
    /** Sparse dwdp temporary storage (dimension: `ndwdp`) */
    mutable SUNMatrixWrapper dwdp_;
    
    /** Dense Mass matrix (dimension: `nx_solver` x `nx_solver`) */
    mutable SUNMatrixWrapper M_;
    
    /**
     * Temporary storage of `dxdotdp_full` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     *  type `CSC_MAT`)
     */
    mutable SUNMatrixWrapper dxdotdp_full;

    /**
     * Temporary storage of `dxdotdp_explicit` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz:  `ndxdotdp_explicit`,
     *  type `CSC_MAT`)
     */
    mutable SUNMatrixWrapper dxdotdp_explicit;

    /**
     * Temporary storage of `dxdotdp_implicit` data across functions,
     * Python-only
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     * type `CSC_MAT`)
     */
    mutable SUNMatrixWrapper dxdotdp_implicit;
    
    /**
     * Temporary storage of `dxdotdx_explicit` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz: 'nxdotdotdx_explicit',
     *  type `CSC_MAT`)
     */
    mutable SUNMatrixWrapper dxdotdx_explicit;

    /**
     * Temporary storage of `dxdotdx_implicit` data across functions,
     * Python-only
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     * type `CSC_MAT`)
     */
    mutable SUNMatrixWrapper dxdotdx_implicit;

    /**
     * Temporary storage of `dxdotdp` data across functions, Matlab only
     * (dimension: `nplist` x `nx_solver`, row-major)
     */
    AmiVectorArray dxdotdp {0, 0};

    /** Current observable (dimension: `nytrue`) */
    mutable std::vector<realtype> my_;

    /** Current event measurement (dimension: `nztrue`) */
    mutable std::vector<realtype> mz_;

    /** Sparse observable derivative of data likelihood, only used if
     * `pythonGenerated` == `true` (dimension `nytrue`, `nJ` x `ny`, row-major)
     */
    mutable std::vector<SUNMatrixWrapper> dJydy_;

    /** Observable derivative of data likelihood, only used if
     * `pythonGenerated` == `false` (dimension `nJ` x `ny` x `nytrue`,
     * row-major)
     */
    mutable std::vector<realtype> dJydy_matlab_;

    /** Observable sigma derivative of data likelihood
     * (dimension nJ x ny x nytrue, row-major)
     */
    mutable std::vector<realtype> dJydsigma_;

    /** State derivative of data likelihood
     * (dimension nJ x nx_solver, row-major)
     */
    mutable std::vector<realtype> dJydx_;

    /** Parameter derivative of data likelihood for current timepoint
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

    /** temporary storage of w data across functions (dimension: nw) */
    mutable std::vector<realtype> w_;

    /** temporary storage for flattened sx,
     * (dimension: nx_solver x nplist, row-major)
     */
    mutable std::vector<realtype> sx_;

    /** temporary storage for x_rdata (dimension: nx_rdata) */
    mutable std::vector<realtype> x_rdata_;

    /** temporary storage for sx_rdata slice (dimension: nx_rdata) */
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

    /** original user-provided, possibly scaled parameter array (dimension: np)
     */
    std::vector<realtype> original_parameters_;

    /** index indicating to which event an event output belongs */
    std::vector<int> z2event_;

    /** state initialization (size nx_solver) */
    std::vector<realtype> x0data_;

    /** sensitivity initialization (size nx_rdata x nplist, row-major) */
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

    /** flag indicating whether steadystate sensitivities are to be computed
     *  via FSA when steadyStateSimulation is used
     */
    SteadyStateSensitivityMode steadystate_sensitivity_mode_ {SteadyStateSensitivityMode::newtonOnly};

    /** flag indicating whether reinitialization of states depending on
     *  fixed parameters is activated
     */
    bool reinitialize_fixed_parameter_initial_states_ {false};

    /** Indicates whether the result of every call to `Model::f*` should be
     * checked for finiteness */
    bool always_check_finite_ {false};

  private:
    /** Sparse dwdp implicit temporary storage (dimension: `ndwdp`) */
    mutable std::vector<SUNMatrixWrapper> dwdp_hierarchical_;

    /** Sparse dwdw temporary storage (dimension: `ndwdw`) */
    mutable SUNMatrixWrapper dwdw_;
    
    /** Sparse dwdx implicit temporary storage (dimension: `ndwdx`) */
    mutable std::vector<SUNMatrixWrapper> dwdx_hierarchical_;
    
    /** Recursion */
    int w_recursion_depth_ {0};
};

bool operator==(const Model &a, const Model &b);

} // namespace amici

#endif // AMICI_MODEL_H
