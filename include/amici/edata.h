#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include "amici/defines.h"

#include <vector>

namespace amici {

class Model;
class ReturnData;

/** @brief ExpData carries all information about experimental or condition-specific data */
class ExpData {

  public:
    /** default constructor */
    ExpData();

    /**
     * constructor that only initializes dimensions
     * @param nytrue
     * @param nztrue
     * @param nmaxevent
     */
    ExpData(int nytrue, int nztrue, int nmaxevent);

    /**
     * constructor that initializes timepoints from vectors
     *
     * @param nytrue               (dimension: scalar)
     * @param nztrue               (dimension: scalar)
     * @param nmaxevent            (dimension: scalar)
     * @param ts                   (dimension: nt)
     */
    ExpData(int nytrue, int nztrue, int nmaxevent,
            std::vector<realtype> ts);
    
    
    /**
     * constructor that initializes timepoints and data from vectors
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
    ExpData(int nytrue, int nztrue, int nmaxevent,
            std::vector<realtype> ts,
            std::vector<realtype> observedData,
            std::vector<realtype> observedDataStdDev,
            std::vector<realtype> observedEvents,
            std::vector<realtype> observedEventsStdDev);

    /**
     * constructor that initializes with Model
     *
     * @param model pointer to model specification object
     */
    ExpData(const Model &model);
    
    /**
     * constructor that initializes with returnData, adds
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y scalar standard deviations for all observables
     * @param sigma_z scalar standard deviations for all event observables
     */
    ExpData(const ReturnData &rdata, realtype sigma_y, realtype sigma_z);
    
    /**
     * constructor that initializes with returnData, adds
     *
     * @param rdata return data pointer with stored simulation results
     * @param sigma_y vector of standard deviations for observables (dimension: nytrue or nt x nytrue, row-major)
     * @param sigma_z vector of standard deviations for event observables (dimension: nztrue or nmaxevent x nztrue, row-major)
     */
    ExpData(const ReturnData &rdata, std::vector<realtype> sigma_y, std::vector<realtype> sigma_z);

    /**
     * @brief Copy constructor
     * @param other object to copy from
     */
    ExpData (const ExpData &other);
    

    ~ExpData() = default;
    
    /**
     * number of timepoints
     *
     * @return number of timepoints
     */
     const int nt() const;

    /**
     * set function that copies data from input to ExpData::ts
     *
     * @param ts timepoints
     */
    void setTimepoints(const std::vector<realtype> &ts);
    
    /**
     * get function that copies data from ExpData::ts to output
     *
     * @return ExpData::ts
     */
    std::vector<realtype> const& getTimepoints() const;
    
    /**
     * get function that returns timepoint at index
     *
     * @param it timepoint index
     * @return timepoint timepoint at index
     */
    realtype getTimepoint(int it) const;
    
    /**
     * set function that copies data from input to ExpData::my
     *
     * @param observedData observed data (dimension: nt x nytrue, row-major)
     */
    void setObservedData(const std::vector<realtype> &observedData);
    
    /**
     * set function that copies observed data for specific observable
     *
     * @param observedData observed data (dimension: nt)
     * @param iy oberved data index
     */
    void setObservedData(const std::vector<realtype> &observedData, int iy);
    
    /**
     * get function that checks whether data at specified indices has been set
     *
     * @param it time index
     * @param iy observable index
     * @return boolean specifying if data was set
     */
    bool isSetObservedData(int it, int iy) const;
    
    /**
     * get function that copies data from ExpData::observedData to output
     *
     * @return observed data (dimension: nt x nytrue, row-major)
     */
    std::vector<realtype> const& getObservedData() const;
    
    /**
     * get function that returns a pointer to observed data at index
     *
     * @param it timepoint index
     * @return pointer to observed data at index (dimension: nytrue)
     */
    const realtype *getObservedDataPtr(int it) const;
    
    /**
     * set function that copies data from input to ExpData::observedDataStdDev
     *
     * @param observedDataStdDev standard deviation of observed data (dimension: nt x nytrue, row-major)
     */
    void setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev);
    
    /**
     * set function that sets all ExpData::observedDataStdDev to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     */
    void setObservedDataStdDev(const realtype stdDev);
    
    /**
     * set function that copies standard deviation of observed data for specific observable
     *
     * @param observedDataStdDev standard deviation of observed data (dimension: nt)
     * @param iy observed data index
     */
    void setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev, int iy);
    
    /**
     * set function that sets all standard deviation of a specific observable to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     * @param iy observed data index
     */
    void setObservedDataStdDev(const realtype stdDev, int iy);
    
    /**
     * get function that checks whether standard deviation of data at specified indices has been set
     *
     * @param it time index
     * @param iy observable index
     * @return boolean specifying if standard deviation of data was set
     */
    bool isSetObservedDataStdDev(int it, int iy) const;
    
    /**
     * get function that copies data from ExpData::observedDataStdDev to output
     *
     * @return standard deviation of observed data
     */
    std::vector<realtype> const& getObservedDataStdDev() const;
    
    /**
     * get function that returns a pointer to standard deviation of observed data at index
     *
     * @param it timepoint index
     * @return pointer to standard deviation of observed data at index
     */
    const realtype *getObservedDataStdDevPtr(int it) const;
    
    /**
     * set function that copies observed event data from input to ExpData::observedEvents
     *
     * @param observedEvents observed data (dimension: nmaxevent x nztrue, row-major)
     */
    void setObservedEvents(const std::vector<realtype> &observedEvents);
    
    /**
     * set function that copies observed event data for specific event observable
     *
     * @param observedEvents observed data (dimension: nmaxevent)
     * @param iz observed event data index
     */
    void setObservedEvents(const std::vector<realtype> &observedEvents, int iz);
    
    /**
     * get function that checks whether event data at specified indices has been set
     *
     * @param ie event index
     * @param iz event observable index
     * @return boolean specifying if data was set
     */
    bool isSetObservedEvents(int ie, int iz) const;
    
    /**
     * get function that copies data from ExpData::mz to output
     *
     * @return observed event data
     */
    std::vector<realtype> const& getObservedEvents() const;
    
    /**
     * get function that returns a pointer to observed data at ieth occurence
     *
     * @param ie event occurence
     * @return pointer to observed event data at ieth occurence
     */
    const realtype *getObservedEventsPtr(int ie) const;
    
    /**
     * set function that copies data from input to ExpData::observedEventsStdDev
     *
     * @param observedEventsStdDev standard deviation of observed event data
     */
    void setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev);
    
    /**
     * set function that sets all ExpData::observedDataStdDev to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     */
    void setObservedEventsStdDev(const realtype stdDev);
    
    /**
     * set function that copies standard deviation of observed data for specific observable
     *
     * @param observedEventsStdDev standard deviation of observed data (dimension: nmaxevent)
     * @param iz observed data index
     */
    void setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz);
    
    /**
     * set function that sets all standard deviation of a specific observable to the input value
     *
     * @param stdDev standard deviation (dimension: scalar)
     * @param iz observed data index
     */
    void setObservedEventsStdDev(const realtype stdDev, int iz);
    
    /**
     * get function that checks whether standard deviation of even data at specified indices has been set
     *
     * @param ie event index
     * @param iz event observable index
     * @return boolean specifying if standard deviation of event data was set
     */
    bool isSetObservedEventsStdDev(int ie, int iz) const;
    
    /**
     * get function that copies data from ExpData::observedEventsStdDev to output
     *
     * @return standard deviation of observed event data
     */
    std::vector<realtype> const& getObservedEventsStdDev() const;
    
    /**
     * get function that returns a pointer to standard deviation of observed event data at ieth occurence
     *
     * @param ie event occurence
     * @return pointer to standard deviation of observed event data at ieth occurence
     */
    const realtype *getObservedEventsStdDevPtr(int ie) const;
    
    /** number of observables */
    const int nytrue;
    /** number of event observables */
    const int nztrue;
    /** maximal number of event occurences */
    const int nmaxevent;
    
    /** condition-specific parameters of size Model::nk() or empty */
    std::vector<realtype> fixedParameters;
    /** condition-specific parameters for pre-equilibration of size Model::nk() or empty.
     * Overrides Solver::newton_preeq */
    std::vector<realtype> fixedParametersPreequilibration;
    /** condition-specific parameters for pre-simulation of size Model::nk() or empty. */
    std::vector<realtype> fixedParametersPresimulation;
    /**
     * @brief duration of pre-simulation
     * if this is > 0, presimualation will be performed from (model->t0 - t_presim) to model->t0
     * using the fixedParameters in fixedParametersPresimulation
     */
    realtype t_presim = 0;

protected:
    
    /**
     * checker for dimensions of input observedData or observedDataStdDev
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void checkDataDimension(std::vector<realtype> input, const char *fieldname) const;
    
    /**
     * checker for dimensions of input observedEvents or observedEventsStdDev
     *
     * @param input vector input to be checked
     * @param fieldname name of the input
     */
    void checkEventsDimension(std::vector<realtype> input, const char *fieldname) const;
    
    /**
     * checks input vector of sigmas for not strictly positive values
     *
     * @param sigmaVector vector input to be checked
     * @param vectorName name of the input
     */
    void checkSigmaPositivity(std::vector<realtype> sigmaVector, const char *vectorName) const;
    
    /**
     * checks input scalar sigma for not strictly positive value
     *
     * @param sigma input to be checked
     * @param sigmaName name of the input
     */
    void checkSigmaPositivity(realtype sigma, const char *sigmaName) const;
    
    /** observation timepoints (dimension: nt) */
    std::vector<realtype> ts;
    
    /** observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> observedData;
    /** standard deviation of observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> observedDataStdDev;
    
    /** observed events (dimension: nmaxevents x nztrue, row-major) */
    std::vector<realtype> observedEvents;
    /** standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)*/
    std::vector<realtype> observedEventsStdDev;
};

} // namespace amici

#endif /* AMICI_EDATA_H */
