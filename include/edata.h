
#define my edata->am_my
#define ysigma edata->am_ysigma
#define mz edata->am_mz
#define zsigma edata->am_zsigma

#ifndef _MY_EDATA
#define _MY_EDATA

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

/** @brief struct that carries all information about experimental data */
typedef struct {
    /** observed data */
    double *am_my; 
    /** standard deviation of observed data */
    double *am_ysigma; 
    
    /** observed events */
    double *am_mz; 
    /** standard deviation of observed events */
    double *am_zsigma; 
    
	} *ExpData;
#endif /* _MY_EDATA */
