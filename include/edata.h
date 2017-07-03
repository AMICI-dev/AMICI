#include "include/udata.h"
#ifndef _MY_EDATA
#define _MY_EDATA

/** @brief struct that carries all information about experimental data */
class ExpData {

public:
    /**
     * @brief Default constructor
     */
    ExpData();
    ExpData(const UserData *udata);
    ~ExpData();
    
    void setDefaults();
    
    /** observed data */
    double *my;
    /** standard deviation of observed data */
    double *sigmay;
    
    /** observed events */
    double *mz;
    /** standard deviation of observed events */
    double *sigmaz;
};

#endif /* _MY_EDATA */
