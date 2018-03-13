#ifndef AMICI_RDATA_H
#define AMICI_RDATA_H

#include "amici/defines.h"

#include <vector>

namespace amici {
class Model;
class ReturnData;
class Solver;
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
 * NOTE: multidimensional arrays are stored in row-major order
 * (FORTRAN-style)
 */
class ReturnData {
  public:
    ReturnData();

    ReturnData(Solver const& solver, const Model *model);

    void invalidate(const realtype t);
    void invalidateLLH();

    void 
    applyChainRuleFactorToSimulationResults(const Model *model);

    ~ReturnData() = default;

    /** timepoints (dimension: nt) */
    std::vector<realtype> ts;

    /** time derivative (dimension: nx) */
    std::vector<realtype> xdot;

    /** Jacobian of differential equation right hand side (dimension: nx x nx,
     * row-major) */
    std::vector<realtype> J;

    /** event output (dimension: nmaxevent x nz, row-major) */
    std::vector<realtype> z;

    /** event output sigma standard deviation (dimension: nmaxevent x nz,
     * row-major) */
    std::vector<realtype> sigmaz;

    /** parameter derivative of event output (dimension: nmaxevent x nz,
     * row-major) */
    std::vector<realtype> sz;

    /** parameter derivative of event output standard deviation (dimension:
     * nmaxevent x nz, row-major)  */
    std::vector<realtype> ssigmaz;

    /** event trigger output (dimension: nmaxevent x nz, row-major)*/
    std::vector<realtype> rz;

    /** parameter derivative of event trigger output (dimension: nmaxevent x nz
     * x nplist, row-major) */
    std::vector<realtype> srz;

    /** second order parameter derivative of event trigger output (dimension:
     * nmaxevent x nztrue x nplist x nplist, row-major) */
    std::vector<realtype> s2rz;

    /** state (dimension: nt x nx, row-major) */
    std::vector<realtype> x;

    /** parameter derivative of state (dimension: nt x nx x nplist,
     * row-major) */
    std::vector<realtype> sx;

    /** observable (dimension: nt x ny, row-major) */
    std::vector<realtype> y;

    /** observable standard deviation (dimension: nt x ny, row-major) */
    std::vector<realtype> sigmay;

    /** parameter derivative of observable (dimension: nt x ny x nplist,
     * row-major) */
    std::vector<realtype> sy;

    /** parameter derivative of observable standard deviation (dimension: nt x
     * ny x nplist, row-major) */
    std::vector<realtype> ssigmay;
    
    /** observable (dimension: nt x ny, row-major) */
    std::vector<realtype> res;

    /** parameter derivative of residual (dimension: nt x ny x nplist,
     * row-major) */
    std::vector<realtype> sres;

    /** number of integration steps forward problem (dimension: nt) */
    std::vector<int> numsteps;

    /** number of integration steps backward problem (dimension: nt) */
    std::vector<int> numstepsB;

    /** number of right hand side evaluations forward problem (dimension: nt) */
    std::vector<int> numrhsevals;

    /** number of right hand side evaluations backwad problem (dimension: nt) */
    std::vector<int> numrhsevalsB;

    /** number of error test failures forward problem (dimension: nt) */
    std::vector<int> numerrtestfails;

    /** number of error test failures backwad problem (dimension: nt) */
    std::vector<int> numerrtestfailsB;

    /** number of linear solver convergence failures forward problem (dimension:
     * nt) */
    std::vector<int> numnonlinsolvconvfails;

    /** number of linear solver convergence failures backwad problem (dimension:
     * nt) */
    std::vector<int> numnonlinsolvconvfailsB;

    /** employed order forward problem (dimension: nt) */
    std::vector<int> order;

    /** flag indicating success of Newton solver */
    int newton_status = 0;

    /** computation time of the Newton solver [s] */
    double newton_time = 0.0;

    /** number of Newton steps for steady state problem */
    std::vector<int> newton_numsteps;

    /** number of linear steps by Newton step for steady state problem */
    std::vector<int> newton_numlinsteps;

    /** preequilibration steady state found be Newton solver */
    std::vector<realtype> x0;

    /** preequilibration sensitivities found be Newton solver */
    std::vector<realtype> sx0;

    /** likelihood value (double[1]) */
    realtype llh = 0.0;

    /** chi2 value (double[1]) */
    realtype chi2 = 0.0;

    /** parameter derivative of likelihood (dimension: nplist) */
    std::vector<realtype> sllh;

    /** second order parameter derivative of likelihood (dimension: (nJ-1) x
     * nplist, row-major) */
    std::vector<realtype> s2llh;

    /** status code (double[1]) */
    int status = 0;

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
    std::vector<AMICI_parameter_scaling> pscale;
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
