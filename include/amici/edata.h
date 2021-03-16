#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/misc.h"
#include "amici/simulation_parameters.h"

#include <vector>

namespace amici {

class Model;
class ReturnData;

/**
 * @brief ExpData carries all information about experimental or
 * condition-specific data
 */
class ExpData : public SimulationParameters {

  public:
    /**
     * @brief default constructor
     */
    ExpData() = default;

    /**
     * @brief Copy constructor, needs to be declared to be generated in
     * swig
     */
    ExpData(const ExpData &) = default;

    /**
     * @brief constructor that only initializes dimensions
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
     * @param fixedParameters Model constants (dimension: nk)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts,
            std::vector<realtype> fixedParameters);

    /**
     * @brief constructor that initializes timepoints and data from vectors
     *
     * @param nytrue Number of observables
     * @param nztrue Number of event outputs
     * @param nmaxevent Maximal number of events to track
     * @param ts Timepoints (dimension: nt)
     * @param observedData observed data (dimension: nt x nytrue, row-major)
     * @param observedDataStdDev standard deviation of observed data
     * (dimension: nt x nytrue, row-major)
     * @param observedEvents observed events
     * (dimension: nmaxevents x nztrue, row-major)
     * @param observedEventsStdDev standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts,
            std::vector<realtype> const &observedData,
            std::vector<realtype> const &observedDataStdDev,
            std::vector<realtype> const &observedEvents,
            std::vector<realtype> const &observedEventsStdDev);

    /**
     * @brief constructor that initializes with Model
     *
     * @param model pointer to model specification object
     */
    explicit ExpData(const Model &model);

    /**
     * @brief constructor that initializes with returnData, adds noise according
     * to specified sigmas
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y scalar standard deviations for all observables
     * @param sigma_z scalar standard deviations for all event observables
     */
    ExpData(const ReturnData &rdata, realtype sigma_y, realtype sigma_z);

    /**
     * @brief constructor that initializes with returnData, adds noise according
     * to specified sigmas
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y vector of standard deviations for observables
     * (dimension: nytrue or nt x nytrue, row-major)
     * @param sigma_z vector of standard deviations for event observables
     * (dimension: nztrue or nmaxevent x nztrue, row-major)
     */
    ExpData(const ReturnData &rdata, std::vector<realtype> sigma_y,
            std::vector<realtype> sigma_z);

    ~ExpData() = default;

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
     * @brief Set function that copies data from input to ExpData::ts
     *
     * @param ts timepoints
     */
    void setTimepoints(const std::vector<realtype> &ts);

    /**
     * @brief get function that copies data from ExpData::ts to output
     *
     * @return ExpData::ts
     */
    std::vector<realtype> const &getTimepoints() const;

    /**
     * @brief get function that returns timepoint at index
     *
     * @param it timepoint index
     *
     * @return timepoint timepoint at index
     */
    realtype getTimepoint(int it) const;

    /**
     * @brief set function that copies data from input to ExpData::my
     *
     * @param observedData observed data (dimension: nt x nytrue, row-major)
     */
    void setObservedData(const std::vector<realtype> &observedData);

    /**
     * @brief set function that copies observed data for specific observable
     *
     * @param observedData observed data (dimension: nt)
     * @param iy observed data index
     */
    void setObservedData(const std::vector<realtype> &observedData, int iy);

    /**
     * @brief get function that checks whether data at specified indices has
     * been set
     *
     * @param it time index
     * @param iy observable index
     *
     * @return boolean specifying if data was set
     */
    bool isSetObservedData(int it, int iy) const;

    /**
     * @brief get function that copies data from ExpData::observedData to output
     *
     * @return observed data (dimension: nt x nytrue, row-major)
     */
    std::vector<realtype> const &getObservedData() const;

    /**
     * @brief get function that returns a pointer to observed data at index
     *
     * @param it timepoint index
     *
     * @return pointer to observed data at index (dimension: nytrue)
     */
    const realtype *getObservedDataPtr(int it) const;

    /**
     * @brief set function that copies data from input to
     * ExpData::observedDataStdDev
     *
     * @param observedDataStdDev standard deviation of observed data (dimension:
     * nt x nytrue, row-major)
     */
    void setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev);

    /**
     * @brief set function that sets all ExpData::observedDataStdDev to the
     * input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     */
    void setObservedDataStdDev(realtype stdDev);

    /**
     * @brief set function that copies standard deviation of observed data for
     * specific observable
     *
     * @param observedDataStdDev standard deviation of observed data (dimension:
     * nt)
     * @param iy observed data index
     */
    void setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev,
                               int iy);

    /**
     * @brief set function that sets all standard deviation of a specific
     * observable to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     * @param iy observed data index
     */
    void setObservedDataStdDev(realtype stdDev, int iy);

    /**
     * @brief get function that checks whether standard deviation of data at
     * specified indices has been set
     *
     * @param it time index
     * @param iy observable index
     * @return boolean specifying if standard deviation of data was set
     */
    bool isSetObservedDataStdDev(int it, int iy) const;

    /**
     * @brief get function that copies data from ExpData::observedDataStdDev to
     * output
     *
     * @return standard deviation of observed data
     */
    std::vector<realtype> const &getObservedDataStdDev() const;

    /**
     * @brief get function that returns a pointer to standard deviation of
     * observed data at index
     *
     * @param it timepoint index
     * @return pointer to standard deviation of observed data at index
     */
    const realtype *getObservedDataStdDevPtr(int it) const;

    /**
     * @brief set function that copies observed event data from input to
     * ExpData::observedEvents
     *
     * @param observedEvents observed data (dimension: nmaxevent x nztrue,
     * row-major)
     */
    void setObservedEvents(const std::vector<realtype> &observedEvents);

    /**
     * @brief set function that copies observed event data for specific event
     * observable
     *
     * @param observedEvents observed data (dimension: nmaxevent)
     * @param iz observed event data index
     */
    void setObservedEvents(const std::vector<realtype> &observedEvents, int iz);

    /**
     * @brief get function that checks whether event data at specified indices
     * has been set
     *
     * @param ie event index
     * @param iz event observable index
     * @return boolean specifying if data was set
     */
    bool isSetObservedEvents(int ie, int iz) const;

    /**
     * @brief get function that copies data from ExpData::mz to output
     *
     * @return observed event data
     */
    std::vector<realtype> const &getObservedEvents() const;

    /**
     * @brief get function that returns a pointer to observed data at ieth
     * occurrence
     *
     * @param ie event occurrence
     *
     * @return pointer to observed event data at ieth occurrence
     */
    const realtype *getObservedEventsPtr(int ie) const;

    /**
     * @brief set function that copies data from input to
     * ExpData::observedEventsStdDev
     *
     * @param observedEventsStdDev standard deviation of observed event data
     */
    void
    setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev);

    /**
     * @brief set function that sets all ExpData::observedDataStdDev to the
     * input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     */
    void setObservedEventsStdDev(realtype stdDev);

    /**
     * @brief set function that copies standard deviation of observed data for
     * specific observable
     *
     * @param observedEventsStdDev standard deviation of observed data
     * (dimension: nmaxevent)
     * @param iz observed data index
     */
    void
    setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev,
                            int iz);

    /**
     * @brief set function that sets all standard deviation of a specific
     * observable to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     * @param iz observed data index
     */
    void setObservedEventsStdDev(realtype stdDev, int iz);

    /**
     * @brief get function that checks whether standard deviation of even data
     * at specified indices has been set
     *
     * @param ie event index
     * @param iz event observable index
     * @return boolean specifying if standard deviation of event data was set
     */
    bool isSetObservedEventsStdDev(int ie, int iz) const;

    /**
     * @brief get function that copies data from ExpData::observedEventsStdDev
     * to output
     *
     * @return standard deviation of observed event data
     */
    std::vector<realtype> const &getObservedEventsStdDev() const;

    /**
     * @brief get function that returns a pointer to standard deviation of
     * observed event data at ie-th occurrence
     *
     * @param ie event occurrence
     *
     * @return pointer to standard deviation of observed event data at ie-th
     * occurrence
     */
    const realtype *getObservedEventsStdDevPtr(int ie) const;

  protected:
    /**
     * @brief resizes observedData, observedDataStdDev, observedEvents and
     * observedEventsStdDev
     */
    void applyDimensions();

    /**
     * @brief resizes observedData and observedDataStdDev
     */
    void applyDataDimension();

    /**
     * @brief resizes observedEvents and observedEventsStdDev
     */
    void applyEventDimension();

    /**
     * @brief checker for dimensions of input observedData or observedDataStdDev
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void checkDataDimension(std::vector<realtype> const &input,
                            const char *fieldname) const;

    /**
     * @brief checker for dimensions of input observedEvents or
     * observedEventsStdDev
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void checkEventsDimension(std::vector<realtype> const &input,
                              const char *fieldname) const;

    /** @brief number of observables */
    int nytrue_{0};

    /** @brief number of event observables */
    int nztrue_{0};

    /** @brief maximal number of event occurrences */
    int nmaxevent_{0};

    /** @brief observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> observed_data_;

    /**
     * @brief standard deviation of observed data (dimension: nt x nytrue,
     * row-major)
     */
    std::vector<realtype> observed_data_std_dev_;

    /**
     * @brief observed events (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> observed_events_;

    /**
     * @brief standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> observed_events_std_dev_;
};

/**
 * @brief checks input vector of sigmas for not strictly positive values
 *
 * @param sigmaVector vector input to be checked
 * @param vectorName name of the input
 */
void checkSigmaPositivity(std::vector<realtype> const &sigmaVector,
                          const char *vectorName);

/**
 * @brief checks input scalar sigma for not strictly positive value
 *
 * @param sigma input to be checked
 * @param sigmaName name of the input
 */
void checkSigmaPositivity(realtype sigma, const char *sigmaName);

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
        Model *model, const ExpData *edata = nullptr,
        FixedParameterContext fpc = FixedParameterContext::simulation);

    ConditionContext &operator=(const ConditionContext &other) = delete;

    ~ConditionContext();

    /**
     * @brief Apply condition-specific settings from edata to the
     * constructor-supplied model, not changing the settings which were
     * backed-up in the constructor call.
     *
     * @param edata
     * @param fpc flag indicating which fixedParameter from edata to apply
     */
    void applyCondition(const ExpData *edata,
                        FixedParameterContext fpc);

    /**
     * @brief Restore original settings on constructor-supplied amici::Model.
     * Will be called during destruction. Explicit call is generally not
     * necessary.
     */
    void restore();

  private:
    Model *model_ = nullptr;
    std::vector<realtype> original_x0_;
    std::vector<realtype> original_sx0_;
    std::vector<realtype> original_parameters_;
    std::vector<realtype> original_fixed_parameters_;
    std::vector<realtype> original_timepoints_;
    std::vector<int> original_parameter_list_;
    std::vector<amici::ParameterScaling> original_scaling_;
    bool original_reinitialize_fixed_parameter_initial_states_;
    std::vector<int> original_reinitialization_state_idxs;
};

} // namespace amici

#endif /* AMICI_EDATA_H */
