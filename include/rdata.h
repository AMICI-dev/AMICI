#define tsdata rdata->am_tsdata
#define xdotdata rdata->am_xdotdata
#define dxdotdpdata rdata->am_dxdotdpdata
#define dydxdata rdata->am_dydxdata
#define dydpdata rdata->am_dydpdata
#define Jdata rdata->am_Jdata
#define rootdata rdata->am_rootdata
#define rootSdata rdata->am_rootSdata
#define rootS2data rdata->am_rootS2data
#define rootvaldata rdata->am_rootvaldata
#define rootvalSdata rdata->am_rootvalSdata
#define rootvalS2data rdata->am_rootvalS2data
#define xdata rdata->am_xdata
#define xSdata rdata->am_xSdata
#define ydata rdata->am_ydata
#define ySdata rdata->am_ySdata

#define numstepsdata rdata->am_numstepsdata
#define numrhsevalsdata rdata->am_numrhsevalsdata
#define orderdata rdata->am_orderdata

#define numstepsSdata rdata->am_numstepsSdata
#define numrhsevalsSdata rdata->am_numrhsevalsSdata

#define llhdata rdata->am_llhdata
#define llhSdata rdata->am_llhSdata
#define llhS2data rdata->am_llhS2data
#define chi2data rdata->am_chi2data


#ifndef _MY_RDATA
#define _MY_RDATA

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

/** @brief struct that stores all data which is later returned by the mex function */
typedef struct {

    /** timepoints */
    double *am_tsdata; 
    /** time derivative */
    double *am_xdotdata; 
    /** parameter derivative of time derivative */
    double *am_dxdotdpdata; 
    /** state derivative of observables */
    double *am_dydxdata; 
    /** parameter derivative of observables */
    double *am_dydpdata; 
    /** Jacobian of differential equation right hand side */
    double *am_Jdata; 
    /** events */
    double *am_rootdata;
    /** parameter derivative of events */
    double *am_rootSdata;
    /** second order parameter derivative of events */
    double *am_rootS2data; 
    /** value of event function */
    double *am_rootvaldata;
    /** parameter derivative of event function */
    double *am_rootvalSdata;
    /** second order parameter derivative of event function */
    double *am_rootvalS2data; /** returned vector */
    /** state */
    double *am_xdata;
    /** parameter derivative of state */
    double *am_xSdata; 
    /** observable */
    double *am_ydata;
    /** parameter derivative of observable */
    double *am_ySdata;
    
    /** number of integration steps forward problem */
    double *am_numstepsdata;
    /** number of integration steps backward problem */
    double *am_numstepsSdata; 
    /** number of right hand side evaluations forward problem */
    double *am_numrhsevalsdata;
    /** number of right hand side evaluations backwad problem */
    double *am_numrhsevalsSdata; 
    /** employed order forward problem */
    double *am_orderdata;
    
    /** likelihood value */
    double *am_llhdata;
    /** chi2 value */
    double *am_chi2data;
    /** parameter derivative of likelihood */
    double *am_llhSdata;
    /** second order parameter derivative of likelihood */
    double *am_llhS2data; 
    
	} *ReturnData;
#endif /* _MY_RDATA */
