#ifndef _MY_EDATA
#define _MY_EDATA

/** @brief struct that carries all information about experimental data */
typedef struct edata {
    /** observed data */
    double *am_my; 
    /** standard deviation of observed data */
    double *am_ysigma; 
    
    /** observed events */
    double *am_mz; 
    /** standard deviation of observed events */
    double *am_zsigma; 
    
	} ExpData;

#ifdef AMICI_WITHOUT_MATLAB
void freeExpData(ExpData *edata);
#endif

#endif /* _MY_EDATA */
