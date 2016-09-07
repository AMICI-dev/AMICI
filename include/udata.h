#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

#define qpositivex udata->am_qpositivex

#define plist udata->am_plist
#define np udata->am_np
#define ny udata->am_ny
#define nytrue udata->am_nytrue
#define nx udata->am_nx
#define nxtrue udata->am_nxtrue
#define nz udata->am_nz
#define nztrue udata->am_nztrue
#define ne udata->am_ne
#define nt udata->am_nt
#define nw udata->am_nw
#define ndwdx udata->am_ndwdx
#define ndwdp udata->am_ndwdp
#define nnz udata->am_nnz
#define nmaxevent udata->am_nmaxevent
#define ng udata->am_ng

#define p udata->am_p
#define k udata->am_k

#define tstart udata->am_tstart
#define ts udata->am_ts

#define pbar udata->am_pbar
#define xbar udata->am_xbar

#define idlist udata->am_idlist

#define sensi udata->am_sensi
#define atol udata->am_atol
#define rtol udata->am_rtol
#define maxsteps udata->am_maxsteps

#define ism udata->am_ism
#define sensi_meth udata->am_sensi_meth

#define linsol udata->am_linsol
#define interpType udata->am_interpType

#define lmm udata->am_lmm
#define iter udata->am_iter

#define stldet udata->am_stldet

#define ubw udata->am_ubw
#define lbw udata->am_lbw

#define b_x0 udata->am_bx0
#define x0data udata->am_x0data
#define b_sx0 udata->am_bsx0
#define sx0data udata->am_sx0data

#define ordering udata->am_ordering

#define tmp_J udata->am_J
#define tmp_dxdotdp udata->am_dxdotdp
#define w_tmp udata->am_w
#define dwdx_tmp udata->am_dwdx
#define dwdp_tmp udata->am_dwdp
#define M_tmp udata->am_M
#define dfdx_tmp udata->am_dfdx
#define stau_tmp udata->am_stau

#define z2event udata->am_z2event
#define h udata->am_h

#ifndef _MY_UDATA
#define _MY_UDATA

/** @brief struct that stores all user provided data */
typedef struct {
    /** positivity flag */
    double *am_qpositivex;
    
    /** parameter reordering */
    int    *am_plist;
    /** number of parameters */
    int    am_np;
    /** number of observables */
    int    am_ny;
    /** number of observables in the unaugmented system */
    int    am_nytrue;
    /** number of states */
    int    am_nx;
    /** number of states in the unaugmented system */
    int    am_nxtrue;
    /** number of event outputs */
    int    am_nz;
    /** number of event outputs in the unaugmented system */
    int    am_nztrue;
    /** number of events */
    int    am_ne;
    /** number of timepoints */
    int    am_nt;
    /** dimension of the augmented objective function for 2nd order ASA */
    int    am_ng;    
    /** number of common expressions */
    int    am_nw;
    /** number of derivatives of common expressions wrt x */
    int    am_ndwdx;
    /** number of derivatives of common expressions wrt p */
    int    am_ndwdp;
    /** number of nonzero entries in jacobian */
    int    am_nnz;
    /** maximal number of events to track */
    int    am_nmaxevent;
    
    /** parameter array */
    double *am_p;
    /** constants array */
    double *am_k;
    
    /** starting time */
    double am_tstart;
    /** timepoints */
    double *am_ts;
    
    /** scaling of parameters */
    double *am_pbar;
    /** scaling of states */
    double *am_xbar;
    
    /** flag array for DAE equations */
    double *am_idlist;
    
    /** flag indicating whether sensitivities are supposed to be computed */
    int am_sensi;
    /** absolute tolerances for integration */
    double am_atol;
    /** relative tolerances for integration */
    double am_rtol;
    /** maximum number of allowed integration steps */
    int am_maxsteps;
    
    /** internal sensitivity method */
    /*!
     * a flag used to select the sensitivity solution method. Its value can be CV SIMULTANEOUS or CV STAGGERED. Only applies for Forward Sensitivities.
     */
    int am_ism;
    
    /** method for sensitivity computation */
    /*!
     * CW_FSA for forward sensitivity analysis, CW_ASA for adjoint sensitivity analysis
     */
    int am_sensi_meth;
    /** linear solver specification */
    int am_linsol;
    /** interpolation type */
    /*!
     * specifies the interpolation type for the forward problem solution which is then used for the backwards problem. can be either CV_POLYNOMIAL or CV_HERMITE
     */
    int am_interpType;
    
    /** linear multistep method */
    /*!
     * specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.
     */
    int am_lmm;
    
    /** nonlinear solver */
    /*!
     * specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.
     */
    int am_iter;
    
    /** flag controlling stability limit detection */
    booleantype am_stldet;
    
    /** upper bandwith of the jacobian */
    int am_ubw;
    /** lower bandwith of the jacobian */
    int am_lbw;
    
    /** flag for state initialisation */
    /*!
     * flag which determines whether analytic state initialisation or provided initialisation should be used
     */
    booleantype am_bx0;
    
    /** flag for sensitivity initialisation */
    /*!
     * flag which determines whether analytic sensitivities initialisation or provided initialisation should be used
     */
    booleantype am_bsx0;
    
    /** state initialisation */
    double *am_x0data;
    
    /** sensitivity initialisation */
    double *am_sx0data;
    
    
    /** state ordering */
    int am_ordering;
    
    /** index indicating to which event an event output belongs */
    double *am_z2event;
    
    /** flag indicating whether a certain heaviside function should be active or not */
    double *am_h;
    
    /** tempory storage of Jacobian data across functions */
    SlsMat am_J;
    /** tempory storage of dxdotdp data across functions */
    realtype *am_dxdotdp;
    /** tempory storage of w data across functions */
    realtype *am_w;
    /** tempory storage of dwdx data across functions */
    realtype *am_dwdx;
    /** tempory storage of dwdp data across functions */
    realtype *am_dwdp;
    /** tempory storage of M data across functions */
    realtype *am_M;
    /** tempory storage of dfdx data across functions */
    realtype *am_dfdx;
    /** tempory storage of stau data across functions */
    realtype *am_stau;

    /** flag indicating whether a NaN in dxdotdp has been reported */
    booleantype am_nan_dxdotdp;
    /** flag indicating whether a NaN in J has been reported */
    booleantype am_nan_J;
    /** flag indicating whether a NaN in JSparse has been reported */
    booleantype am_nan_JSparse;
    /** flag indicating whether a NaN in xdot has been reported */
    booleantype am_nan_xdot;
    /** flag indicating whether a NaN in xBdot has been reported */
    booleantype am_nan_xBdot;
    /** flag indicating whether a NaN in qBdot has been reported */
    booleantype am_nan_qBdot;
    
    
} *UserData;
#endif /* _MY_UDATA */
