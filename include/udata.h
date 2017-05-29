#ifndef _MY_UDATA
#define _MY_UDATA

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

typedef enum AMI_parameter_scaling_TAG {
    AMI_SCALING_NONE, AMI_SCALING_LN, AMI_SCALING_LOG10
} AMI_parameter_scaling;

typedef enum AMI_o2mode_TAG {
    AMI_O2MODE_NONE, AMI_O2MODE_FULL, AMI_O2MODE_DIR
} AMI_o2mode;

typedef enum AMI_sensi_order_TAG {
    AMI_SENSI_ORDER_NONE, AMI_SENSI_ORDER_FIRST, AMI_SENSI_ORDER_SECOND
} AMI_sensi_order;

typedef enum AMI_sensi_meth_TAG {
    AMI_SENSI_NONE, AMI_SENSI_FSA, AMI_SENSI_ASA, AMI_SENSI_SS
} AMI_sensi_meth;


/** @brief struct that stores all user provided data */
class UserData {

public:
    UserData(int np,
             int nx, int nxtrue,
             int nk,
             int ny, int nytrue,
             int nz, int nztrue,
             int ne, int ng,
             int nw, int ndwdx, int ndwdp, int nnz,
             int ubw, int lbw,
             AMI_parameter_scaling pscale,
             AMI_o2mode o2mode
             );

    virtual ~UserData();

    /* Model dimensions */
    /** total number of model parameters */
    const int    np;
    /** number of fixed parameters */
    const int    nk;
    /** number of observables */
    const int    ny;
    /** number of observables in the unaugmented system */
    const int    nytrue;
    /** number of states */
    const int    nx;
    /** number of states in the unaugmented system */
    const int    nxtrue;
    /** number of event outputs */
    const int    nz;
    /** number of event outputs in the unaugmented system */
    const int    nztrue;
    /** number of events */
    const int    ne;
    /** number of common expressions */
    const int    nw;
    /** number of derivatives of common expressions wrt x */
    const int    ndwdx;
    /** number of derivatives of common expressions wrt p */
    const int    ndwdp;
    /** number of nonzero entries in jacobian */
    const int    nnz;
    /** flag indicating whether for sensi == 2 directional or full second order derivative will be computed */
    const AMI_o2mode o2mode;
    /** dimension of the augmented objective function for 2nd order ASA */
    const int    ng;
    /** upper bandwith of the jacobian */
    const int ubw;
    /** lower bandwith of the jacobian */
    const int lbw;

    /* Options */

    /** maximal number of events to track */
    int    nmaxevent;

    /** positivity flag */
    double *qpositivex;

    /** parameter selection and reordering */
    int    *plist;
    /** number of parameters in plist */
    int    nplist;

    /** number of timepoints */
    int    nt;

    /** parametrization of parameters p */
    AMI_parameter_scaling pscale;

    /** parameter array */
    double *p;
    /** constants array */
    double *k;
    
    /** starting time */
    double tstart;
    /** timepoints */
    double *ts;
    
    /** scaling of parameters */
    double *pbar;
    /** scaling of states */
    double *xbar;
    
    /** flag array for DAE equations */
    double *idlist;
    
    /** flag indicating whether sensitivities are supposed to be computed */
    AMI_sensi_order sensi;
    /** absolute tolerances for integration */
    double atol;
    /** relative tolerances for integration */
    double rtol;
    /** maximum number of allowed integration steps */
    int maxsteps;
    
    /** internal sensitivity method */
    /*!
     * a flag used to select the sensitivity solution method. Its value can be CV SIMULTANEOUS or CV STAGGERED. Only applies for Forward Sensitivities.
     */
    int ism;
    
    /** method for sensitivity computation */
    /*!
     * CW_FSA for forward sensitivity analysis, CW_ASA for adjoint sensitivity analysis
     */
    AMI_sensi_meth sensi_meth;
    /** linear solver specification */
    int linsol;
    /** interpolation type */
    /*!
     * specifies the interpolation type for the forward problem solution which is then used for the backwards problem. can be either CV_POLYNOMIAL or CV_HERMITE
     */
    int interpType;
    
    /** linear multistep method */
    /*!
     * specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.
     */
    int lmm;
    
    /** nonlinear solver */
    /*!
     * specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.
     */
    int iter;
    
    /** flag controlling stability limit detection */
    booleantype stldet;
    
    

    /** state initialisation */
    double *x0data;
    
    /** sensitivity initialisation */
    double *sx0data;
    
    
    /** state ordering */
    int ordering;
    
    /** index indicating to which event an event output belongs */
    double *z2event;
    
    /** flag indicating whether a certain heaviside function should be active or not */
    double *h;
    
    /** tempory storage of Jacobian data across functions */
    SlsMat J;
    /** tempory storage of dxdotdp data across functions */
    realtype *dxdotdp;
    /** tempory storage of w data across functions */
    realtype *w;
    /** tempory storage of dwdx data across functions */
    realtype *dwdx;
    /** tempory storage of dwdp data across functions */
    realtype *dwdp;
    /** tempory storage of M data across functions */
    realtype *M;
    /** tempory storage of dfdx data across functions */
    realtype *dfdx;
    /** tempory storage of stau data across functions */
    realtype *stau;

    /** flag indicating whether a NaN in dxdotdp has been reported */
    booleantype nan_dxdotdp;
    /** flag indicating whether a NaN in J has been reported */
    booleantype nan_J;
    /** flag indicating whether a NaN in JSparse has been reported */
    booleantype nan_JSparse;
    /** flag indicating whether a NaN in xdot has been reported */
    booleantype nan_xdot;
    /** flag indicating whether a NaN in xBdot has been reported */
    booleantype nan_xBdot;
    /** flag indicating whether a NaN in qBdot has been reported */
    booleantype nan_qBdot;
        
};

#ifdef AMICI_WITHOUT_MATLAB
void printUserData(UserData *udata);
#endif

#endif /* _MY_UDATA */
