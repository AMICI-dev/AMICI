#define tsdata rdata->am_tsdata
#define xdotdata rdata->am_xdotdata
#define dxdotdpdata rdata->am_dxdotdpdata
#define dydxdata rdata->am_dydxdata
#define dydpdata rdata->am_dydpdata
#define Jdata rdata->am_Jdata
#define zdata rdata->am_zdata
#define sigmazdata rdata->am_sigmazdata
#define szdata rdata->am_szdata
#define ssigmazdata rdata->am_ssigmazdata
#define rzdata rdata->am_rzdata
#define srzdata rdata->am_srzdata
#define s2rzdata rdata->am_s2rzdata
#define xdata rdata->am_xdata
#define sxdata rdata->am_sxdata
#define ydata rdata->am_ydata
#define sigmaydata rdata->am_sigmaydata
#define sydata rdata->am_sydata
#define ssigmaydata rdata->am_ssigmaydata

#define numstepsdata rdata->am_numstepsdata
#define numrhsevalsdata rdata->am_numrhsevalsdata
#define orderdata rdata->am_orderdata

#define numstepsSdata rdata->am_numstepsSdata
#define numrhsevalsSdata rdata->am_numrhsevalsSdata

#define llhdata rdata->am_llhdata
#define sllhdata rdata->am_sllhdata
#define s2llhdata rdata->am_s2llhdata
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
    /** event output */
    double *am_zdata;
    /** event output sigma standard deviation */
    double *am_sigmazdata;
    /** parameter derivative of event output */
    double *am_szdata;
    /** parameter derivative of event output standard deviation */
    double *am_ssigmazdata;
    /** event trigger output */
    double *am_rzdata;
    /** parameter derivative of event trigger output */
    double *am_srzdata;
    /** second order parameter derivative of event trigger output */
    double *am_s2rzdata;
    /** state */
    double *am_xdata;
    /** parameter derivative of state */
    double *am_sxdata;
    /** observable */
    double *am_ydata;
    /** observable standard deviation */
    double *am_sigmaydata;
    /** parameter derivative of observable */
    double *am_sydata;
    /** parameter derivative of observable standard deviation */
    double *am_ssigmaydata;
    
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
    double *am_sllhdata;
    /** second order parameter derivative of likelihood */
    double *am_s2llhdata;
    
	} *ReturnData;

#ifdef AMICI_WITHOUT_MATLAB
void freeReturnData(ReturnData rdata);
#endif

#endif /* _MY_RDATA */
