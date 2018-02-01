#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include <vector>
#include <include/amici_defines.h>

namespace amici {

class Model;

/** @brief struct that carries all information about experimental data */
class ExpData {

  public:
    /** default constructor */
    ExpData();
    ExpData(Model *model);

    /**
     * @brief ExpData is currently not copyable
     * @param other object to copy from
     */
    ExpData (const ExpData &other) = delete;

    void setObservedData(const double *observedData);
    void setObservedDataStdDev(const double *observedDataStdDev);
    void setObservedEvents(const double *observedEvents);
    void setObservedEventsStdDev(const double *observedEventsStdDev);

    ~ExpData();

    /** observed data (dimension: nt x nytrue, column-major) */
    std::vector<realtype> my;
    /** standard deviation of observed data (dimension: nt x nytrue, column-major) */
    std::vector<realtype> sigmay;

    /** observed events (dimension: nmaxevents x nztrue, column-major) */
    std::vector<realtype> mz;
    /** standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, column-major)*/
    std::vector<realtype> sigmaz;
    
    /** number of observables */
    const int nytrue;
    /** number of event observables */
    const int nztrue;
    /** number of timepoints */
    const int nt;
    /** maximal number of event occurences */
    const int nmaxevent;
};

} // namespace amici

#endif /* _MY_EDATA */
