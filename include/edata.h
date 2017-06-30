#include "include/udata.h"
#ifndef _MY_EDATA
#define _MY_EDATA

#ifdef AMICI_WITHOUT_MATLAB
#define mxArray double
#else
#include <mex.h>
#endif

/** @brief struct that carries all information about experimental data */
class ExpData {

public:
    /**
     * @brief Default constructor
     */
    ExpData(const UserData *udata);
    ~ExpData();

    int expDataFromMatlabCall(const mxArray *prhs[], const UserData *udata);
    int expDataFromCppCall(const char* hdffile, const UserData *udata);
    
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
