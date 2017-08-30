#ifndef _MY_UDATA
#define _MY_UDATA
#include "include/amici_defines.h"

#include <cmath>

/** @brief struct that stores all user provided data
 * NOTE: multidimensional arrays are expected to be stored in column-major order
 * (FORTRAN-style)
 */
class UserData {

  public:
    /**
     * @brief UserData
     * @param np total number of model parameters
     * @param nk number of f ixed parameters
     * @param nx number of states
     */
    UserData(int np, int nk, int nx);

    /**
     * @brief Default constructor for testing and serialization
     */
    UserData();

    int unscaleParameters(double *bufferUnscaled) const;

    /**
     * @brief setTimepoints
     * @param timepoints
     * @param numTimepoints
     */
    void setTimepoints(const double *timepoints, int numTimepoints);

    /**
     * @brief setParameters
     * @param parameters
     */
    void setParameters(const double *parameters);

    /**
     * @brief setConstants
     * @param constants
     */
    void setConstants(const double *constants);

    /**
     * @brief setPlist set parameter selection and ordering
     * @param plist
     */
    void setPlist(const double *plist, int nplist);
    void setPlist(const int *plist, int nplist);

    /**
     * @brief setPbar. Must not be called before setPlist
     * @param parameterScaling
     */
    void setPbar(const double *parameterScaling);

    /**
     * @brief setStateInitialization
     * @param stateInitialization
     */
    void setStateInitialization(const double *stateInitialization);

    /**
     * @brief setSensitivityInitialization
     * @param sensitivityInitialization
     */
    void setSensitivityInitialization(const double *sensitivityInitialization);

    ~UserData();

    /* Options */

    /** maximal number of events to track */
    int nmaxevent = 10;

    /** positivity flag (size nx) */
    double *qpositivex = nullptr;

    /** parameter selection and reordering (size nplist) */
    int *plist = nullptr;

    /** number of parameters in plist */
    int nplist = 0;

    /** number of timepoints */
    int nt = 0;

    /** parameter array (size np) */
    double *p = nullptr;

    /** constants array (size nk) */
    double *k = nullptr;

    /** parameter transformation of p */
    AMICI_parameter_scaling pscale = AMICI_SCALING_NONE;

    /** starting time */
    double tstart = 0.0;

    /** timepoints (size nt) */
    double *ts = nullptr;

    /** scaling of parameters (size nplist) */
    double *pbar = nullptr;

    /** scaling of states (size nx) */
    double *xbar = nullptr;

    /** flag indicating whether sensitivities are supposed to be computed */
    AMICI_sensi_order sensi = AMICI_SENSI_ORDER_NONE;

    /** absolute tolerances for integration */
    double atol = 1e-16;

    /** relative tolerances for integration */
    double rtol = 1e-8;

    /** maximum number of allowed integration steps */
    int maxsteps = 0;

    /** maximum number of allowed Newton steps for steady state computation */
    int newton_maxsteps = 0;

    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    int newton_maxlinsteps = 0;

    /** Preequilibration of model via Newton solver? */
    int newton_preeq = false;

    /** Which preconditioner is to be used in the case of iterative linear
     * Newton solvers */
    int newton_precon = 1;

    /** internal sensitivity method flag used to select the sensitivity solution
     * method. Its value can be CV SIMULTANEOUS or CV STAGGERED. Only applies
     * for Forward Sensitivities. */
    int ism = 1;

    /** method for sensitivity computation */
    AMICI_sensi_meth sensi_meth = AMICI_SENSI_FSA;

    /** linear solver specification */
    int linsol = 9;

    /** interpolation type for the forward problem solution which
     * is then used for the backwards problem. can be either CV_POLYNOMIAL or
     * CV_HERMITE
     */
    int interpType = 1;

    /** specifies the linear multistep method and may be one of two possible
     * values: CV ADAMS or CV BDF.
     */
    int lmm = 2;

    /**
     * specifies the type of nonlinear solver iteration and may be either CV
     * NEWTON or CV FUNCTIONAL.
     */
    int iter = 2;

    /** flag controlling stability limit detection */
    booleantype stldet = true;

    /** state initialisation (size np) */
    double *x0data = nullptr;

    /** sensitivity initialisation (size nx * nplist) */
    double *sx0data = nullptr;

    /** state ordering */
    int ordering = 0;

    /** function to print the contents of the UserData object */
    void print() const;

    /** total number of model parameters */
    const int np;
    /** number of fixed parameters */
    const int nk;
    /** number of states */
    const int nx;
};

#endif /* _MY_UDATA */
