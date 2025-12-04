#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include "amici/defines.h"
#include "amici/misc.h"
#include "amici/simulation_parameters.h"

#include <string>
#include <vector>

namespace amici {

class Model;
class ReturnData;

/**
 * @brief ExpData carries all information about experimental or
 * condition-specific data.
 */
class ExpData : public SimulationParameters {

  public:
    /**
     * @brief Default constructor.
     */
    ExpData() = default;

    /**
     * @brief Copy constructor.
     */
    // needs to be declared to be wrapped by SWIG
    ExpData(ExpData const&) = default;

    /**
     * @brief Constructor that only initializes dimensions.
     *
     * @param nytrue Number of observables
     * @param nztrue Number of event outputs
     * @param nmaxevent Maximal number of events to track
     */
    ExpData(int nytrue, int nztrue, int nmaxevent);

    /**
     * @brief constructor that initializes timepoints from vectors
     *
     * @param nytrue Number of observables
     * @param nztrue Number of event outputs
     * @param nmaxevent Maximal number of events to track
     * @param ts Timepoints (dimension: nt)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts);

    /**
     * @brief constructor that initializes timepoints and fixed parameters from
     * vectors
     *
     * @param nytrue Number of observables
     * @param nztrue Number of event outputs
     * @param nmaxevent Maximal number of events to track
     * @param ts Timepoints (dimension: nt)
     * @param fixed_parameters Model variables excluded from sensitivity
     * analysis (dimension: nk)
     */
    ExpData(
        int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts,
        std::vector<realtype> fixed_parameters
    );

    /**
     * @brief constructor that initializes timepoints and data from vectors
     *
     * @param nytrue Number of observables
     * @param nztrue Number of event outputs
     * @param nmaxevent Maximal number of events to track
     * @param ts Timepoints (dimension: nt)
     * @param my measurements (dimension: nt x nytrue, row-major)
     * @param sigma_y noise scale of measurements
     * (dimension: nt x nytrue, row-major)
     * @param mz event measurements
     * (dimension: nmaxevents x nztrue, row-major)
     * @param sigma_z noise scale of event measurements
     * (dimension: nmaxevents x nztrue, row-major)
     */
    ExpData(
        int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts,
        std::vector<realtype> const& my,
        std::vector<realtype> const& sigma_y,
        std::vector<realtype> const& mz,
        std::vector<realtype> const& sigma_z
    );

    /**
     * @brief constructor that initializes with Model
     *
     * @param model pointer to model specification object
     */
    explicit ExpData(Model const& model);

    /**
     * @brief Constructor that initializes with ReturnData, adds normally
     * distributed noise according to specified sigmas.
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y scalar noise scales for all observables
     * @param sigma_z scalar noise scales for all event observables
     * @param seed Seed for the random number generator. If a negative number
     * is passed, a random seed is used.
     */
    ExpData(
        ReturnData const& rdata, realtype sigma_y, realtype sigma_z,
        int seed = -1
    );

    /**
     * @brief Constructor that initializes with ReturnData, adds normally
     * distributed noise according to specified sigmas.
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y vector of noise scales for observables
     * (dimension: nytrue or nt x nytrue, row-major)
     * @param sigma_z vector of noise scales for event observables
     * (dimension: nztrue or nmaxevent x nztrue, row-major)
     * @param seed Seed for the random number generator. If a negative number
     * is passed, a random seed is used.
     */
    ExpData(
        ReturnData const& rdata, std::vector<realtype> sigma_y,
        std::vector<realtype> sigma_z, int seed = -1
    );

    ~ExpData() = default;

    friend inline bool operator==(ExpData const& lhs, ExpData const& rhs);

    /**
     * @brief number of observables of the non-augmented model
     *
     * @return number of observables of the non-augmented model
     */
    int nytrue() const;

    /**
     * @brief number of event observables of the non-augmented model
     *
     * @return number of event observables of the non-augmented model
     */
    int nztrue() const;

    /**
     * @brief maximal number of events to track
     *
     * @return maximal number of events to track
     */
    int nmaxevent() const;

    /**
     * @brief number of timepoints
     *
     * @return number of timepoints
     */
    int nt() const;

    /**
     * @brief Set output ts.
     *
     * If the number of timepoint increases, this will grow the
     * observation/sigma matrices and fill new entries with NaN.
     * If the number of ts decreases, this will shrink the
     * observation/sigma matrices.
     *
     * Note that the mapping from ts to measurements will not be
     * preserved. E.g., say there are measurements at t = 2, and this
     * function is called with [1, 2], then the old measurements will belong to
     * t = 1.
     *
     * @param ts ts
     */
    void set_timepoints(std::vector<realtype> const& ts);

    /**
     * @brief Get output timepoints.
     *
     * @return ExpData::ts
     */
    std::vector<realtype> const& get_timepoints() const;

    /**
     * @brief Get timepoint for a specific index.
     *
     * @param it timepoint index
     *
     * @return timepoint timepoint at index
     */
    realtype get_timepoint(int it) const;

    /**
     * @brief Set all measurements.
     *
     * @param my measurements (dimension: nt x nytrue, row-major)
     */
    void set_measurements(std::vector<realtype> const& my);

    /**
     * @brief Set measurements for a specific observable.
     *
     * @param my measurements (dimension: nt)
     * @param iy observable index
     */
    void set_measurements(std::vector<realtype> const& my, int iy);

    /**
     * @brief Check whether a measurement is defined at the given indices.
     *
     * @param it time index
     * @param iy observable index
     *
     * @return true if a value is set for the specified indices; otherwise false
     */
    bool is_set_measurement(int it, int iy) const;

    /**
     * @brief Get all measurements.
     *
     * @return measurements (dimension: nt x nytrue, row-major)
     */
    std::vector<realtype> const& get_measurements() const;

    /**
     * @brief Get measurements for a specific timepoint.
     *
     * @param it timepoint index
     *
     * @return pointer to measurements at index (dimension: nytrue)
     */
    realtype const* get_measurements_ptr(int it) const;

    /**
     * @brief Set noise scales for all measurements.
     *
     * @param sigma noise scales (dimension: nt x nytrue, row-major)
     */
    void set_noise_scales(std::vector<realtype> const& sigma);

    /**
     * @brief Set identical noise scales for all measurements.
     *
     * @param sigma noise scale (dimension: scalar)
     */
    void set_noise_scales(realtype sigma);

    /**
     * @brief Set measurement noise scales of for a observable.
     *
     * @param sigma noise scales (dimension: nt)
     * @param iy observable index
     */
    void set_noise_scales(std::vector<realtype> const& sigma, int iy);

    /**
     * @brief Set all noise scales for a specific observable to the
     * input value.
     *
     * @param sigma noise scale (dimension: scalar)
     * @param iy observable index
     */
    void set_noise_scales(realtype sigma, int iy);

    /**
     * @brief Check whether a noise scale is defined at the given indices.
     *
     * @param it time index
     * @param iy observable index
     * @return true if a value is set for the specified indices; otherwise false
     */
    bool is_set_noise_scale(int it, int iy) const;

    /**
     * @brief Get measurement noise scales.
     *
     * @return noise scales of measurements
     */
    std::vector<realtype> const& get_noise_scales() const;

    /**
     * @brief Get pointer to measurement noise scales.
     *
     * @param it timepoint index
     * @return pointer to noise scales of measurements at index
     */
    realtype const* get_noise_scales_ptr(int it) const;

    /**
     * @brief Set event measurements.
     *
     * @param mz event measurements (dimension: nmaxevent x nztrue,
     * row-major)
     */
    void set_event_measurements(std::vector<realtype> const& mz);

    /**
     * @brief Set event measurements for a specific event observable.
     *
     * @param mz event measurements (dimension: nmaxevent)
     * @param iz event observable index
     */
    void
    set_event_measurements(std::vector<realtype> const& mz, int iz);

    /**
     * @brief Check whether an event measurement is defined at the given indices.
     *
     * @param ie event index
     * @param iz event observable index
     * @return true if a value is set for the specified indices; otherwise false
     */
    bool is_set_event_measurement(int ie, int iz) const;

    /**
     * @brief Get all event measurements.
     *
     * @return event measurements
     */
    std::vector<realtype> const& get_event_measurements() const;

    /**
     * @brief Get pointer to event measurements at ie-th occurrence.
     *
     * @param ie event occurrence
     *
     * @return pointer to event measurements at ie-th occurrence
     */
    realtype const* get_event_measurements_ptr(int ie) const;

    /**
     * @brief Set noise scales of event measurements.
     *
     * @param sigma noise scales of event measurements
     */
    void set_event_noise_scales(std::vector<realtype> const& sigma);

    /**
     * @brief Set noise scales of all event measurements.
     *
     * @param sigma noise scale (dimension: scalar)
     */
    void set_event_noise_scales(realtype sigma);

    /**
     * @brief Set noise scales for a specific event observable.
     *
     * @param sigma noise scales of observed data
     * (dimension: nmaxevent)
     * @param iz event observable index
     */
    void set_event_noise_scales(std::vector<realtype> const& sigma, int iz);

    /**
     * @brief Set all noise scales for a specific event observable.
     *
     * @param sigma noise scale (dimension: scalar)
     * @param iz event observable index
     */
    void set_event_noise_scales(realtype sigma, int iz);

    /**
     * @brief Check whether an event noise scale is defined at the given indices.
     *
     * @param ie event occurence
     * @param iz event observable index
     * @return true if a value is set for the specified indices; otherwise false
     */
    bool is_set_event_noise_scale(int ie, int iz) const;

    /**
     * @brief Get noise scale of observed event data.
     *
     * @return noise scale of observed event data
     */
    std::vector<realtype> const& get_event_noise_scales() const;

    /**
     * @brief Get pointer to noise scale of
     * observed event data at ie-th occurrence.
     *
     * @param ie event occurrence
     *
     * @return pointer to noise scale of observed event data at ie-th
     * occurrence
     */
    realtype const* get_event_noise_scales_ptr(int ie) const;

    /**
     * @brief Set all observations and their noise scales to NaN.
     *
     * Useful, e.g., after calling ExpData::setTimepoints.
     */
    void clear_observations();

    /**
     * @brief Arbitrary (not necessarily unique) identifier.
     */
    std::string id;

  protected:
    /**
     * @brief resizes measurements_, noise_scales_, event_measurements_ and
     * event_noise_scales_
     */
    void apply_dimensions();

    /**
     * @brief resizes measurements_ and noise_scales_
     */
    void apply_data_dimension();

    /**
     * @brief resizes event_measurements_ and event_noise_scales_
     */
    void apply_event_dimension();

    /**
     * @brief checker for dimensions of input measurements_ or noise_scales_
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void check_data_dimension(
        std::vector<realtype> const& input, char const* fieldname
    ) const;

    /**
     * @brief checker for dimensions of input event_measurements_ or
     * event_noise_scales_
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void check_events_dimension(
        std::vector<realtype> const& input, char const* fieldname
    ) const;

    /** @brief number of observables */
    int nytrue_{0};

    /** @brief number of event observables */
    int nztrue_{0};

    /** @brief maximal number of event occurrences */
    int nmaxevent_{0};

    /** @brief observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> measurements_;

    /**
     * @brief noise scale of observed data (dimension: nt x nytrue,
     * row-major)
     */
    std::vector<realtype> noise_scales_;

    /**
     * @brief observed events (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> event_measurements_;

    /**
     * @brief noise scale of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> event_noise_scales_;
};

/**
 * @brief Equality operator
 *
 * @param lhs some object
 * @param rhs another object
 * @return `true`, if both arguments are equal; `false` otherwise.
 */
inline bool operator==(ExpData const& lhs, ExpData const& rhs) {
    return *dynamic_cast<SimulationParameters const*>(&lhs)
               == *dynamic_cast<SimulationParameters const*>(&rhs)
           && lhs.id == rhs.id && lhs.nytrue_ == rhs.nytrue_
           && lhs.nztrue_ == rhs.nztrue_ && lhs.nmaxevent_ == rhs.nmaxevent_
           && is_equal(lhs.measurements_, rhs.measurements_)
           && is_equal(lhs.noise_scales_, rhs.noise_scales_)
           && is_equal(lhs.event_measurements_, rhs.event_measurements_)
           && is_equal(
               lhs.event_noise_scales_, rhs.event_noise_scales_
           );
}

/**
 * @brief checks input vector of sigmas for not strictly positive values
 *
 * @param sigmaVector vector input to be checked
 * @param vectorName name of the input
 */
void check_sigma_positivity(
    std::vector<realtype> const& sigmaVector, char const* vectorName
);

/**
 * @brief checks input scalar sigma for not strictly positive value
 *
 * @param sigma input to be checked
 * @param sigmaName name of the input
 */
void check_sigma_positivity(realtype sigma, char const* sigmaName);

/**
 * @brief The ConditionContext class applies condition-specific amici::Model
 * settings and restores them when going out of scope
 */
class ConditionContext : public ContextManager {
  public:
    /**
     * @brief Apply condition-specific settings from edata to model while
     * keeping a backup of the original values.
     *
     * @param model
     * @param edata
     * @param fpc flag indicating which fixedParameter from edata to apply
     */
    explicit ConditionContext(
        Model* model, ExpData const* edata = nullptr,
        FixedParameterContext fpc = FixedParameterContext::simulation
    );

    ConditionContext& operator=(ConditionContext const& other) = delete;

    ~ConditionContext();

    /**
     * @brief Apply condition-specific settings from edata to the
     * constructor-supplied model, not changing the settings which were
     * backed-up in the constructor call.
     *
     * @param edata
     * @param fpc flag indicating which fixedParameter from edata to apply
     */
    void apply_condition(ExpData const* edata, FixedParameterContext fpc);

    /**
     * @brief Restore original settings on constructor-supplied amici::Model.
     *
     * Will be called during destruction. Explicit call is generally not
     * necessary.
     */
    void restore();

  private:
    Model* model_ = nullptr;
    std::vector<realtype> original_x0_;
    std::vector<realtype> original_sx0_;
    std::vector<realtype> original_parameters_;
    std::vector<realtype> original_fixed_parameters_;
    realtype original_tstart_;
    realtype original_tstart_preeq_;
    std::vector<realtype> original_timepoints_;
    std::vector<int> original_parameter_list_;
    std::vector<ParameterScaling> original_scaling_;
    bool original_reinitialize_fixed_parameter_initial_states_;
    std::vector<int> original_reinitialization_state_idxs;
};

} // namespace amici

#endif /* AMICI_EDATA_H */
