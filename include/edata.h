#ifndef _MY_EDATA
#define _MY_EDATA

class UserData;
class Model;

/** @brief struct that carries all information about experimental data */
class ExpData {

public:
    /**
     * @brief Default constructor
     */
    ExpData();
    /** initialization with UserData and model */
    ExpData(const UserData *udata, Model *model);
    ~ExpData();
    
    /** initialization with default values */
    void setDefaults();

    /** observed data */
    double *my;
    /** standard deviation of observed data */
    double *sigmay;

    /** observed events */
    double *mz;
    /** observed roots */
    double *mrz;
    /** standard deviation of observed events/roots */
    double *sigmaz;
};

#endif /* _MY_EDATA */
