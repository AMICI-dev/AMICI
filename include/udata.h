#ifndef _MY_UDATA
#define _MY_UDATA

#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>

typedef enum AMICI_parameter_scaling_TAG {
    AMICI_SCALING_NONE, AMICI_SCALING_LN, AMICI_SCALING_LOG10
} AMICI_parameter_scaling;

typedef enum AMICI_o2mode_TAG {
    AMICI_O2MODE_NONE, AMICI_O2MODE_FULL, AMICI_O2MODE_DIR
} AMICI_o2mode;

typedef enum AMICI_sensi_order_TAG {
    AMICI_SENSI_ORDER_NONE, AMICI_SENSI_ORDER_FIRST, AMICI_SENSI_ORDER_SECOND
} AMICI_sensi_order;

typedef enum AMICI_sensi_meth_TAG {
    AMICI_SENSI_NONE, AMICI_SENSI_FSA, AMICI_SENSI_ASA, AMICI_SENSI_SS
} AMICI_sensi_meth;


/** @brief struct that stores all user provided data */
class UserData {

public:
    /**
     * @brief Default constructor for testing and serialization
     */
    UserData();

    UserData(int np,
             int nx, int nxtrue,
             int nk,
             int ny, int nytrue,
             int nz, int nztrue,
             int ne, int nJ,
             int nw, int ndwdx, int ndwdp, int nnz,
             int ubw, int lbw,
             AMICI_parameter_scaling pscale,
             AMICI_o2mode o2mode
             );

    void initTemporaryFields();
    void freeTemporaryFields();

    virtual ~UserData();

    /* Model dimensions */
    /** total number of model parameters */
    const int    np;
    /** number of fixed parameters */
    const int    nk;
    /** number of states */
    const int    nx;
    /** number of states in the unaugmented system */
    const int    nxtrue;
    /** number of observables */
    const int    ny;
    /** number of observables in the unaugmented system */
    const int    nytrue;
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
    /** dimension of the augmented objective function for 2nd order ASA */
    const int    nJ;
    /** upper bandwith of the jacobian */
    const int ubw;
    /** lower bandwith of the jacobian */
    const int lbw;
    /** flag indicating whether for sensi == 2 directional or full second order derivative will be computed */
    const AMICI_o2mode o2mode;

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
    AMICI_parameter_scaling pscale;

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
    AMICI_sensi_order sensi;
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
    AMICI_sensi_meth sensi_meth;
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
    realtype *h;
    
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
        
protected:
    void init();
};

#ifdef AMICI_WITHOUT_MATLAB
void printUserData(UserData *udata);
#endif

#endif /* _MY_UDATA */
