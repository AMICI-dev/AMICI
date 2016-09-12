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
    /** old state vector */
    N_Vector am_x_old;
    /** array of state vectors at discontinuities*/
    N_Vector *am_x_disc;
    /** array of differential state vectors at discontinuities*/
    N_Vector *am_xdot_disc;
    /** array of old differential state vectors at discontinuities*/
    N_Vector *am_xdot_old_disc;
    /** differential state vector */
    N_Vector am_dx;
    /** old differential state vector */
    N_Vector am_dx_old;
    /** time derivative state vector */
    N_Vector am_xdot;
    /** old time derivative state vector */
    N_Vector am_xdot_old;
    /** adjoint state vector */
    N_Vector am_xB;
    /** old adjoint state vector */
    N_Vector am_xB_old;
    /** differential adjoint state vector */
    N_Vector am_dxB; 
    /** quadrature state vector */
    N_Vector am_xQB;
    /** old quadrature state vector */
    N_Vector am_xQB_old;
    /** sensitivity state vector array */
    N_Vector *am_sx; 
    /** differential sensitivity state vector array */
    N_Vector *am_sdx; 
    /** index indicating DAE equations vector */
    N_Vector am_id;
    /** Jacobian */
    DlsMat am_Jtmp;
    
    /** parameter derivative of likelihood array */
    realtype *am_llhS0;
    /** data likelihood */
    realtype *am_g;
    /** parameter derivative of data likelihood */
    realtype *am_dgdp;
    /** state derivative of data likelihood */
    realtype *am_dgdx;
    /** event likelihood */
    realtype *am_r;
    /** parameter derivative of event likelihood */
    realtype *am_drdp;
    /** state derivative of event likelihood */
    realtype *am_drdx;
    /** root function likelihood */
    realtype am_rval; 
    /** parameter derivative of root function likelihood */
    realtype *am_drvaldp;
    /** state derivative of root function likelihood */
    realtype *am_drvaldx;
    /** state derivative of event */
    realtype *am_dzdx;
    /** parameter derivative of event */
    realtype *am_dzdp;
    /** parameter derivative of observable */
    realtype *am_dydp;
    /** state derivative of observable */
    realtype *am_dydx;
    /** initial sensitivity of observable */
    realtype *am_yS0;
    /** data standard deviation */
    realtype *am_sigma_y;
    /** parameter derivative of data standard deviation */
    realtype *am_dsigma_ydp;
    /** event standard deviation */
    realtype *am_sigma_z;
    /** parameter derivative of event standard deviation */
    realtype *am_dsigma_zdp;
    
    /** state array */
    realtype *am_x_tmp;
    /** sensitivity state array */
    realtype *am_sx_tmp;
    /** differential state array */
    realtype *am_dx_tmp;
    /** differential sensitivity state array */
    realtype *am_sdx_tmp;
    /** time derivative state array */
    realtype *am_xdot_tmp;
    /** differential adjoint state array */
    realtype *am_xB_tmp;
    /** quadrature state array */
    realtype *am_xQB_tmp;
    /** differential adjoint state array */
    realtype *am_dxB_tmp;
    /** index indicating DAE equations array */
    realtype *am_id_tmp;
    /** temporary storage for heaviside flags to check whether a secondary event has fired */
    realtype *am_h_tmp;

    
    /** array of flags indicating which root has beend found */
    /*!
    array of length nr with the indices of the user functions gi found to have a root. For i = 0, . . . ,nr?1, rootsfound[i]?= 0 if gi has a root, and = 0 if not.
    */
    int *am_rootsfound;
    /** array of index which root has been found */
    int *am_rootidx;
    /** array of number of found roots for a certain event type */
    int *am_nroots;
    /** array of values of the root function */
    double *am_rootvals;
    
    
    /** change in x */
    realtype *am_deltax;
    /** change in sx */
    realtype *am_deltasx;
    /** change in xB */
    realtype *am_deltaxB;
    /** change in qB */
    realtype *am_deltaqB;
 
    
    /** integer for indexing of backwards problems */
    int am_which;
    
    /** array containing the time-points of discontinuities*/
    realtype *am_discs; 
    /** array containing the index of discontinuities */
    realtype *am_irdiscs; 

	} *TempData;
#endif /* _MY_TDATA */
