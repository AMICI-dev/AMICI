#define t tdata->am_t

#define x tdata->am_x
#define dx tdata->am_dx
#define xdot tdata->am_xdot
#define xB tdata->am_xB
#define dxB tdata->am_dxB
#define xQB tdata->am_xQB
#define sx tdata->am_sx
#define sdx tdata->am_sdx
#define Jtmp tdata->am_Jtmp
#define id tdata->am_id

#define llhS0 tdata->am_llhS0
#define g tdata->am_g
#define r tdata->am_r
#define dzdp tdata->am_dzdp
#define dzdx tdata->am_dzdx
#define dgdp tdata->am_dgdp
#define dgdx tdata->am_dgdx
#define drdp tdata->am_drdp
#define drdx tdata->am_drdx
#define drvaldx tdata->am_drvaldx
#define dydp tdata->am_dydp
#define dydx tdata->am_dydx
#define sigma_y tdata->am_sigma_y
#define dsigma_ydp tdata->am_dsigma_ydp
#define sigma_z tdata->am_sigma_z
#define dsigma_zdp tdata->am_dsigma_zdp

#define x_tmp tdata->am_x_tmp
#define dx_tmp tdata->am_dx_tmp
#define xB_tmp tdata->am_xB_tmp
#define dxB_tmp tdata->am_dxB_tmp
#define xQB_tmp tdata->am_xQB_tmp
#define sx_tmp tdata->am_sx_tmp
#define sdx_tmp tdata->am_sdx_tmp
#define xdot_tmp tdata->am_xdot_tmp
#define id_tmp tdata->am_id_tmp

#define rootsfound tdata->am_rootsfound
#define rootvaltmp tdata->am_rootvaltmp
#define rootidx tdata->am_rootidx

#define which tdata->am_which

#define discs tdata->am_discs
#define irdiscs tdata->am_irdiscs


#ifndef _MY_TDATA
#define _MY_TDATA

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

/** @brief struct that provides temporary storage for different variables */
typedef struct {
    /** current time */
    realtype am_t;
    
    
    /** state vector */
    N_Vector am_x; 
    /** differential state vector */
    N_Vector am_dx; 
    /** time derivative state vector */
    N_Vector am_xdot; 
    /** adjoint state vector */
    N_Vector am_xB; 
    /** differential adjoint state vector */
    N_Vector am_dxB; 
    /** quadrature state vector */
    N_Vector am_xQB; 
    /** sensitivity state vector array */
    N_Vector *am_sx; 
    /** differential sensitivity state vector array */
    N_Vector *am_sdx; 
    /** index indicating DAE equations vector */
    N_Vector am_id;
    /** Jacobian */
    DlsMat am_Jtmp;
    
    /** parameter derivative of likelihood array */
    double *am_llhS0;
    /** data likelihood */
    double am_g;
    /** parameter derivative of data likelihood */
    double *am_dgdp;
    /** state derivative of data likelihood */
    double *am_dgdx;
    /** event likelihood */
    double am_r;
    /** parameter derivative of event likelihood */
    double *am_drdp;
    /** state derivative of event likelihood */
    double *am_drdx;
    /** root function likelihood */
    double am_rval; 
    /** parameter derivative of root function likelihood */
    double *am_drvaldp;
    /** state derivative of root function likelihood */
    double *am_drvaldx;
    /** state derivative of event */
    double *am_dzdx;
    /** parameter derivative of event */
    double *am_dzdp;
    /** parameter derivative of observable */
    double *am_dydp;
    /** state derivative of observable */
    double *am_dydx;
    /** initial sensitivity of observable */
    double *am_yS0;
    /** data standard deviation */
    double *am_sigma_y;
    /** parameter derivative of data standard deviation */
    double *am_dsigma_ydp;
    /** event standard deviation */
    double *am_sigma_z;
    /** parameter derivative of event standard deviation */
    double *am_dsigma_zdp;
    
    /** state array */
    double *am_x_tmp;
    /** sensitivity state array */
    double *am_sx_tmp;
    /** differential state array */
    double *am_dx_tmp;
    /** differential sensitivity state array */
    double *am_sdx_tmp;
    /** time derivative state array */
    double *am_xdot_tmp;
    /** differential adjoint state array */
    double *am_xB_tmp;
    /** quadrature state array */
    double *am_xQB_tmp;
    /** differential adjoint state array */
    double *am_dxB_tmp;
    /** index indicating DAE equations array */
    double *am_id_tmp;
    
    /** array of flags indicating which root has beend found */
    /*!
    array of length nr with the indices of the user functions gi found to have a root. For i = 0, . . . ,nr?1, rootsfound[i]?= 0 if gi has a root, and = 0 if not.
    */
    int *am_rootsfound;
    /** array of values of the root function */
    double *am_rootvaltmp; 
    /** array of index which root has been found */
    int *am_rootidx;
 
    
    /** integer for indexing of backwards problems */
    int am_which;
    
    /** array containing the time-points of discontinuities*/
    double *am_discs; 
    /** array containing the index of discontinuities */
    double *am_irdiscs; 

	} *TempData;
#endif /* _MY_TDATA */
