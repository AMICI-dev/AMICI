#ifndef _MY_EDATA
#define _MY_EDATA

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#include <stdbool.h>
#endif

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

EXTERNC void freeExpData(ExpData *edata);

#endif /* _MY_EDATA */
