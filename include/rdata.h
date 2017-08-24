#ifndef _MY_RDATA
#define _MY_RDATA
#include <include/udata.h>
class UserData;
class Model;

/** @brief struct that stores all data which is later returned by the mex
 * function
 *
 * NOTE: MATLAB stores multidimensional arrays in column-major order
 * (FORTRAN-style)
 */
class ReturnData {
  public:
    ReturnData();

    ReturnData(const UserData *udata, const Model *model);

    virtual void setDefaults();

    void invalidate();

    void setLikelihoodSensitivityFirstOrderNaN();
    void setLikelihoodSensitivitySecondOrderNaN();

    int
    applyChainRuleFactorToSimulationResults(const UserData *udata,
                                            const realtype *unscaledParameters);

    virtual ~ReturnData();

    /** timepoints (dimension: nt) */
    double *ts;

    /** time derivative (dimension: nx) */
    double *xdot;

    /** parameter derivative of time derivative (dimension: nx x nplist,
     * column-major) */
    double *dxdotdp;

    /** state derivative of observables (dimension: ny x nx, column-major) */
    double *dydx;

    /** parameter derivative of observables (dimension: ny x nplist,
     * column-major) */
    double *dydp;

    /** Jacobian of differential equation right hand side (dimension: nx x nx,
     * column-major) */
    double *J;

    /** event output (dimension: nmaxevent x nz, column-major) */
    double *z;

    /** event output sigma standard deviation (dimension: nmaxevent x nz,
     * column-major) */
    double *sigmaz;

    /** parameter derivative of event output (dimension: nmaxevent x nz,
     * column-major) */
    double *sz;

    /** parameter derivative of event output standard deviation (dimension:
     * nmaxevent x nz, column-major)  */
    double *ssigmaz;

    /** event trigger output (dimension: nmaxevent x nz, column-major)*/
    double *rz;

    /** parameter derivative of event trigger output (dimension: nmaxevent x nz
     * x nplist, column-major) */
    double *srz;

    /** second order parameter derivative of event trigger output (dimension:
     * nmaxevent x nztrue x nplist x nplist, column-major) */
    double *s2rz;

    /** state (dimension: nt x nx, column-major) */
    double *x;

    /** parameter derivative of state (dimension: nt x nx x nplist,
     * column-major) */
    double *sx;

    /** observable (dimension: nt x ny, column-major) */
    double *y;

    /** observable standard deviation (dimension: nt x ny, column-major) */
    double *sigmay;

    /** parameter derivative of observable (dimension: nt x ny x nplist,
     * column-major) */
    double *sy;

    /** parameter derivative of observable standard deviation (dimension: nt x
     * ny x nplist, column-major) */
    double *ssigmay;

    /** number of integration steps forward problem (dimension: nt) */
    double *numsteps;

    /** number of integration steps backward problem (dimension: nt) */
    double *numstepsB;

    /** number of right hand side evaluations forward problem (dimension: nt) */
    double *numrhsevals;

    /** number of right hand side evaluations backwad problem (dimension: nt) */
    double *numrhsevalsB;

    /** number of error test failures forward problem (dimension: nt) */
    double *numerrtestfails;

    /** number of error test failures backwad problem (dimension: nt) */
    double *numerrtestfailsB;

    /** number of linear solver convergence failures forward problem (dimension:
     * nt) */
    double *numnonlinsolvconvfails;

    /** number of linear solver convergence failures backwad problem (dimension:
     * nt) */
    double *numnonlinsolvconvfailsB;

    /** employed order forward problem (dimension: nt) */
    double *order;

    /** flag indicating success of Newton solver */
    double *newton_status;

    /** computation time of the Newton solver [s] */
    double *newton_time;

    /** number of Newton steps for steady state problem */
    double *newton_numsteps;

    /** number of linear steps by Newton step for steady state problem */
    double *newton_numlinsteps;

    /** steady state found be Newton solver */
    double *xss;

    /** likelihood value (double[1]) */
    double *llh;

    /** chi2 value (double[1]) */
    double *chi2;

    /** parameter derivative of likelihood (dimension: nplist) */
    double *sllh;

    /** second order parameter derivative of likelihood (dimension: (nJ-1) x
     * nplist, column-major) */
    double *s2llh;

    /** status code (double[1]) */
    double *status;

  protected:
    virtual void copyFromUserData(const UserData *udata);

    virtual void initFields();

    virtual void initField1(double **fieldPointer, const char *fieldName,
                            int dim);

    virtual void initField2(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2);

    virtual void initField3(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3);

    virtual void initField4(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3, int dim4);

    /** flag indicating whether memory for fields needs to be freed on destruction */
    bool freeFieldsOnDestruction;

  public:
    /** total number of model parameters */
    const int np;
    /** number of fixed parameters */
    const int nk;
    /** number of states */
    const int nx;
    /** number of states in the unaugmented system */
    const int nxtrue;
    /** number of observables */
    const int ny;
    /** number of observables in the unaugmented system */
    const int nytrue;
    /** number of event outputs */
    const int nz;
    /** number of event outputs in the unaugmented system */
    const int nztrue;
    /** number of events */
    const int ne;
    /** dimension of the augmented objective function for 2nd order ASA */
    const int nJ;

    /** number of parameter for which sensitivities were requested */
    const int nplist;
    /** maximal number of occuring events (for every event type) */
    const int nmaxevent;
    /** number of considered timepoints */
    const int nt;
    /** maximal number of newton iterations for steady state calculation */
    const int newton_maxsteps;
    /** scaling of parameterization (lin,log,log10) */
    const AMICI_parameter_scaling pscale;
    /** flag indicating whether second order sensitivities were requested */
    const AMICI_o2mode o2mode;
    /** sensitivity order */
    const AMICI_sensi_order sensi;
    /** sensitivity method */
    const AMICI_sensi_meth sensi_meth;
};

#endif /* _MY_RDATA */
