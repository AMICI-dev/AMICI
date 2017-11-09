#ifndef AMICI_RDATA_H
#define AMICI_RDATA_H
#include <include/udata.h>

namespace amici {
class Model;
class ReturnData;
}

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::ReturnData &u, const unsigned int version);
}}

namespace amici {

/** @brief class that stores all data which is later returned by the mex
 * function
 *
 * NOTE: multidimensional arrays are stored in column-major order
 * (FORTRAN-style)
 */
class ReturnData {
  public:
    ReturnData();

    ReturnData(const UserData *udata, const Model *model);

    void invalidate(const int t);
    void invalidateLLH();

    void setLikelihoodSensitivityFirstOrderNaN();

    void setLikelihoodSensitivitySecondOrderNaN();

    void 
    applyChainRuleFactorToSimulationResults(const UserData *udata);

    virtual ~ReturnData();

    /** timepoints (dimension: nt) */
    double *ts = nullptr;

    /** time derivative (dimension: nx) */
    double *xdot = nullptr;

    /** Jacobian of differential equation right hand side (dimension: nx x nx,
     * column-major) */
    double *J = nullptr;

    /** event output (dimension: nmaxevent x nz, column-major) */
    double *z = nullptr;

    /** event output sigma standard deviation (dimension: nmaxevent x nz,
     * column-major) */
    double *sigmaz = nullptr;

    /** parameter derivative of event output (dimension: nmaxevent x nz,
     * column-major) */
    double *sz = nullptr;

    /** parameter derivative of event output standard deviation (dimension:
     * nmaxevent x nz, column-major)  */
    double *ssigmaz = nullptr;

    /** event trigger output (dimension: nmaxevent x nz, column-major)*/
    double *rz = nullptr;

    /** parameter derivative of event trigger output (dimension: nmaxevent x nz
     * x nplist, column-major) */
    double *srz = nullptr;

    /** second order parameter derivative of event trigger output (dimension:
     * nmaxevent x nztrue x nplist x nplist, column-major) */
    double *s2rz = nullptr;

    /** state (dimension: nt x nx, column-major) */
    double *x = nullptr;

    /** parameter derivative of state (dimension: nt x nx x nplist,
     * column-major) */
    double *sx = nullptr;

    /** observable (dimension: nt x ny, column-major) */
    double *y = nullptr;

    /** observable standard deviation (dimension: nt x ny, column-major) */
    double *sigmay = nullptr;

    /** parameter derivative of observable (dimension: nt x ny x nplist,
     * column-major) */
    double *sy = nullptr;

    /** parameter derivative of observable standard deviation (dimension: nt x
     * ny x nplist, column-major) */
    double *ssigmay = nullptr;

    /** number of integration steps forward problem (dimension: nt) */
    double *numsteps = nullptr;

    /** number of integration steps backward problem (dimension: nt) */
    double *numstepsB = nullptr;

    /** number of right hand side evaluations forward problem (dimension: nt) */
    double *numrhsevals = nullptr;

    /** number of right hand side evaluations backwad problem (dimension: nt) */
    double *numrhsevalsB = nullptr;

    /** number of error test failures forward problem (dimension: nt) */
    double *numerrtestfails = nullptr;

    /** number of error test failures backwad problem (dimension: nt) */
    double *numerrtestfailsB = nullptr;

    /** number of linear solver convergence failures forward problem (dimension:
     * nt) */
    double *numnonlinsolvconvfails = nullptr;

    /** number of linear solver convergence failures backwad problem (dimension:
     * nt) */
    double *numnonlinsolvconvfailsB = nullptr;

    /** employed order forward problem (dimension: nt) */
    double *order = nullptr;

    /** flag indicating success of Newton solver */
    double *newton_status = nullptr;

    /** computation time of the Newton solver [s] */
    double *newton_time = nullptr;

    /** number of Newton steps for steady state problem */
    double *newton_numsteps = nullptr;

    /** number of linear steps by Newton step for steady state problem */
    double *newton_numlinsteps = nullptr;

    /** preequilibration steady state found be Newton solver */
    double *x0 = nullptr;

    /** preequilibration sensitivities found be Newton solver */
    double *sx0 = nullptr;

    /** likelihood value (double[1]) */
    double *llh = nullptr;

    /** chi2 value (double[1]) */
    double *chi2 = nullptr;

    /** parameter derivative of likelihood (dimension: nplist) */
    double *sllh = nullptr;

    /** second order parameter derivative of likelihood (dimension: (nJ-1) x
     * nplist, column-major) */
    double *s2llh = nullptr;

    /** status code (double[1]) */
    double *status = nullptr;

  protected:

    ReturnData(const UserData *udata, const Model *model,
               bool initializeFields);

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

    /** flag indicating whether memory for fields needs to be freed on
     * destruction */
    bool freeFieldsOnDestruction = true;

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

    /**
     * @brief Serialize ReturnData (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, ReturnData &r, const unsigned int version);
};

} // namespace amici

#endif /* _MY_RDATA */
