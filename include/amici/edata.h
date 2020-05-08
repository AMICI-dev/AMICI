#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/misc.h"

#include <vector>

namespace amici {

class Model;
class ReturnData;

/**
 * @brief ExpData carries all information about experimental or
 * condition-specific data
 */
class ExpData {

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
     * @param nytrue
     * @param nztrue
     * @param nmaxevent
     */
    ExpData(int nytrue, int nztrue, int nmaxevent);

    /**
     * @brief constructor that initializes timepoints from vectors
     *
     * @param nytrue               (dimension: scalar)
     * @param nztrue               (dimension: scalar)
     * @param nmaxevent            (dimension: scalar)
     * @param ts                   (dimension: nt)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts);

    /**
     * @brief constructor that initializes timepoints and fixed parameters from
     * vectors
     *
     * @param nytrue               (dimension: scalar)
     * @param nztrue               (dimension: scalar)
     * @param nmaxevent            (dimension: scalar)
     * @param ts                   (dimension: nt)
     * @param fixedParameters      (dimension: nk)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent, std::vector<realtype> ts,
            std::vector<realtype> fixedParameters);

    /**
     * @brief constructor that initializes timepoints and data from vectors
     *
     * @param nytrue               (dimension: scalar)
     * @param nztrue               (dimension: scalar)
     * @param nmaxevent            (dimension: scalar)
     * @param ts                   (dimension: nt)
     * @param observedData         (dimension: nt x nytrue, row-major)
     * @param observedDataStdDev   (dimension: nt x nytrue, row-major)
     * @param observedEvents       (dimension: nmaxevent x nztrue, row-major)
     * @param observedEventsStdDev (dimension: nmaxevent x nztrue, row-major)
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
     * @brief set function that copies data from input to ExpData::ts
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
     * @param iy oberved data index
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
     * occurence
     *
     * @param ie event occurence
     *
     * @return pointer to observed event data at ieth occurence
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
     * observed event data at ieth occurence
     *
     * @param ie event occurence
     *
     * @return pointer to standard deviation of observed event data at ieth
     * occurence
     */
    const realtype *getObservedEventsStdDevPtr(int ie) const;

    /**
     * @brief condition-specific fixed parameters of size Model::nk() or empty
     */
    std::vector<realtype> fixedParameters;
    /**
     * @brief condition-specific fixed parameters for pre-equilibration of size
     * Model::nk() or empty. Overrides Solver::newton_preeq
     */
    std::vector<realtype> fixedParametersPreequilibration;
    /**
     * @brief condition-specific fixed parameters for pre-simulation of
     * size Model::nk() or empty.
     */
    std::vector<realtype> fixedParametersPresimulation;

    /**
     * @brief condition-specific parameters of size Model::np() or empty
     */
    std::vector<realtype> parameters;
    /**
     * @brief condition-specific initial conditions of size Model::nx() or
     * empty
     */
    std::vector<realtype> x0;
    /**
     * @brief condition-specific initial condition sensitivities of size
     * Model::nx() * Model::nplist(), Model::nx() * ExpDataplist.size(), if
     * ExpData::plist is not empty, or empty
     */
    std::vector<realtype> sx0;
    /**
     * @brief condition-specific parameter scales of size Model::np()
     */
    std::vector<ParameterScaling> pscale;
    /**
     * @brief condition-specific parameter list
     */
    std::vector<int> plist;

    /**
     * @brief duration of pre-simulation
     * if this is > 0, presimualation will be performed from
     * (model->t0 - t_presim) to model->t0 using the fixedParameters in
     * fixedParametersPresimulation
     */
    realtype t_presim = 0;

    /** flag indicating whether reinitialization of states depending on
     *  fixed parameters is activated
     */
    bool reinitializeFixedParameterInitialStates = false;

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
     * @param input vector input to be checkedjupyter_contrib_nbextensions
     * @param fieldname name of the input
     */
    void checkEventsDimension(std::vector<realtype> const &input,
                              const char *fieldname) const;

    /** @brief number of observables */
    int nytrue_{0};

    /** @brief number of event observables */
    int nztrue_{0};

    /** @brief maximal number of event occurences */
    int nmaxevent_{0};

    /** @brief observation timepoints (dimension: nt) */
    std::vector<realtype> ts;

    /** @brief observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> observedData;
    /**
     * @brief standard deviation of observed data (dimension: nt x nytrue,
     * row-major)
     */
    std::vector<realtype> observedDataStdDev;

    /**
     * @brief observed events (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> observedEvents;
    /**
     * @brief standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)
     */
    std::vector<realtype> observedEventsStdDev;
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
     * @param fpc flag indicating which fixedParmeter from edata to apply
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
     * @param fpc flag indicating which fixedParmeter from edata to apply
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
    Model *model = nullptr;
    std::vector<realtype> originalx0;
    std::vector<realtype> originalsx0;
    std::vector<realtype> originalParameters;
    std::vector<realtype> originalFixedParameters;
    std::vector<realtype> originalTimepoints;
    std::vector<int> originalParameterList;
    std::vector<amici::ParameterScaling> originalScaling;
    bool originalReinitializeFixedParameterInitialStates;
};

} // namespace amici

#endif /* AMICI_EDATA_H */
