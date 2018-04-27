#ifndef AMICI_EDATA_H
#define AMICI_EDATA_H

#include "amici/defines.h"

#include <vector>

namespace amici {

class Model;

/** @brief struct that carries all information about experimental data */
class ExpData {

  public:
    /** default constructor */
    ExpData();
    /**
     * @brief ExpData
     * @param nytrue
     * @param nztrue
     * @param nt
     * @param nmaxevent
     */
    ExpData(int nytrue, int nztrue, int nt, int nmaxevent);
    /**
     * @brief ExpData
     * @param nytrue
     * @param nztrue
     * @param nt
     * @param nmaxevent
     * @param my
     * @param sigmay
     * @param mz
     * @param sigmaz
     */
    ExpData(int nytrue, int nztrue, int nt, int nmaxevent,
            std::vector<realtype> const& my,
            std::vector<realtype> const& sigmay,
            std::vector<realtype> const& mz,
            std::vector<realtype> const& sigmaz);
    /**
     * constructor that initializes with Model
     *
     * @param model pointer to model specification object @type Model
     */
    ExpData(const Model &model);

    /**
     * @brief Copy constructor
     * @param other object to copy from
     */
    ExpData (const ExpData &other);

    void setObservedData(const double *observedData);
    void setObservedDataStdDev(const double *observedDataStdDev);
    void setObservedEvents(const double *observedEvents);
    void setObservedEventsStdDev(const double *observedEventsStdDev);

    ~ExpData();

    /** observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> my;
    /** standard deviation of observed data (dimension: nt x nytrue, row-major) */
    std::vector<realtype> sigmay;

    /** observed events (dimension: nmaxevents x nztrue, row-major) */
    std::vector<realtype> mz;
    /** standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, row-major)*/
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
