#ifndef _MY_EDATA
#define _MY_EDATA

class UserData;
class Model;

/** @brief struct that carries all information about experimental data */
class ExpData {

  public:
    /** default constructor */
    ExpData();
    ExpData(const UserData *udata, Model *model);

    void setObservedData(const double *observedData);
    void setObservedDataStdDev(const double *observedDataStdDev);
    void setObservedEvents(const double *observedEvents);
    void setObservedEventsStdDev(const double *observedEventsStdDev);

    ~ExpData();

    /** observed data (dimension: nytrue x nt, column-major) */
    double *my = nullptr;
    /** standard deviation of observed data (dimension: nytrue x nt, column-major) */
    double *sigmay = nullptr;

    /** observed events (dimension: nmaxevents x nztrue, column-major) */
    double *mz = nullptr;
    /** standard deviation of observed events/roots
     * (dimension: nmaxevents x nztrue, column-major)*/
    double *sigmaz = nullptr;
    
    /** number of observables */
    const int nytrue;
    /** number of event observables */
    const int nztrue;
    /** number of timepoints */
    const int nt;
    /** maximal number of event occurences */
    const int nmaxevent;
};

#endif /* _MY_EDATA */
