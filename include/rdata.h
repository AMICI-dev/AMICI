#ifndef _MY_RDATA
#define _MY_RDATA

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/** @brief struct that stores all data which is later returned by the mex function */
typedef struct rdata {

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
    /** steady state */
    double *am_xssdata;
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
    /** flag success of Newton solver */
    double *am_newtdata;
    /** number of Newton steps steady state problem */
    double *am_numstepsNewtdata;
    /** number of linear steps per Newton step steady state problem */
    double *am_numlinstepsNewtdata;
    
    /** likelihood value */
    double *am_llhdata;
    /** chi2 value */
    double *am_chi2data;
    /** parameter derivative of likelihood */
    double *am_sllhdata;
    /** second order parameter derivative of likelihood */
    double *am_s2llhdata;
    
	} ReturnData;

EXTERNC void freeReturnData(ReturnData *rdata);

#endif /* _MY_RDATA */
