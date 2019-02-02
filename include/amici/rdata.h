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

/** @brief Stores all data to be returned by amici::runAmiciSimulation.
 *
 * NOTE: multidimensional arrays are stored in row-major order
 * (FORTRAN-style)
 */
class ReturnData {
  public:
    /**
     * @brief default constructor
     */
    ReturnData();

    /**
     * @brief ReturnData
     * @param ts see amici::Model::ts
     * @param np see amici::Model::np
     * @param nk see amici::Model::nk
     * @param nx see amici::Model::nx_rdata
     * @param nx_solver see amici::Model::nx_solver
     * @param nxtrue see amici::Model::nxtrue_rdata
     * @param ny see amici::Model::ny
     * @param nytrue see amici::Model::nytrue
     * @param nz see amici::Model::nz
     * @param nztrue see amici::Model::nztrue
     * @param ne see amici::Model::ne
     * @param nJ see amici::Model::nJ
     * @param nplist see amici::Model::nplist
     * @param nmaxevent see amici::Model::nmaxevent
     * @param nt see amici::Model::nt
     * @param newton_maxsteps see amici::Solver::newton_maxsteps
     * @param pscale see amici::Model::pscale
     * @param o2mode see amici::Model::o2mode
     * @param sensi see amici::Solver::sensi
     * @param sensi_meth see amici::Solver::sensi_meth
     */
    ReturnData(
            std::vector<realtype> ts,
            int np, int nk, int nx, int nx_solver, int nxtrue, int ny, int nytrue,
            int nz, int nztrue, int ne, int nJ, int nplist, int nmaxevent,
            int nt, int newton_maxsteps, std::vector<ParameterScaling> pscale,
            SecondOrderMode o2mode, SensitivityOrder sensi, SensitivityMethod sensi_meth);

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver
     * @param model pointer to model specification object
     * bool
     */
    ReturnData(Solver const& solver, const Model *model);

    ~ReturnData() = default;

    /**
     * @brief initializeObjectiveFunction
     */
    void initializeObjectiveFunction();

    /**
     * @brief Set likelihood, state variables, outputs and respective
     * sensitivities to NaN (typically after integration failure)
     * @param t time of integration failure
     */
    void invalidate(const realtype t);

    /**
     * @brief Set likelihood and chi2 to NaN
     * (typically after integration failure)
     */
    void invalidateLLH();

    /**
     * @brief Set likelihood sensitivities to NaN
     * (typically after integration failure)
     */
    void invalidateSLLH();

    /**
     * @brief applies the chain rule to account for parameter transformation
     * in the sensitivities of simulation results
     * @param model Model from which the ReturnData was obtained
     */
    void
    applyChainRuleFactorToSimulationResults(const Model *model);

    /** timepoints (dimension: nt) */
    const std::vector<realtype> ts;

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

    /** parameter derivative of state (dimension: nt x nplist x nx,
     * row-major) */
    std::vector<realtype> sx;

    /** observable (dimension: nt x ny, row-major) */
    std::vector<realtype> y;

    /** observable standard deviation (dimension: nt x ny, row-major) */
    std::vector<realtype> sigmay;

    /** parameter derivative of observable (dimension: nt x nplist x ny,
     * row-major) */
    std::vector<realtype> sy;

    /** parameter derivative of observable standard deviation (dimension: nt x
     * nplist x ny, row-major) */
    std::vector<realtype> ssigmay;

    /** observable (dimension: nt*ny, row-major) */
    std::vector<realtype> res;

    /** parameter derivative of residual (dimension: nt*ny x nplist,
     * row-major) */
    std::vector<realtype> sres;

    /** fisher information matrix (dimension: nplist x nplist,
     * row-major) */
    std::vector<realtype> FIM;

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
    double newton_cpu_time = 0.0;

    /** number of Newton steps for steady state problem
     [newton, simulation, newton] (length = 3) */
    std::vector<int> newton_numsteps;

    /** number of linear steps by Newton step for steady state problem. this
     will only be filled for iterative solvers (length = newton_maxsteps * 2) */
    std::vector<int> newton_numlinsteps;

    /** time at which steadystate was reached in the simulation based approach */
    realtype t_steadystate = NAN;

    /** weighted root-mean-square of the rhs when steadystate
     was reached*/
    realtype wrms_steadystate = NAN;

    /** weighted root-mean-square of the rhs when steadystate
     was reached*/
    realtype wrms_sensi_steadystate = NAN;


    /** initial state (dimension: nx) */
    std::vector<realtype> x0;

    /** preequilibration steady state found by Newton solver (dimension: nx) */
    std::vector<realtype> x_ss;

    /** initial sensitivities (dimension: nplist x nx, row-major) */
    std::vector<realtype> sx0;

    /** preequilibration sensitivities found by Newton solver (dimension: nplist x nx, row-major) */
    std::vector<realtype> sx_ss;

    /** loglikelihood value */
    realtype llh = 0.0;

    /** chi2 value */
    realtype chi2 = 0.0;

    /** parameter derivative of loglikelihood (dimension: nplist) */
    std::vector<realtype> sllh;

    /** second order parameter derivative of loglikelihood (dimension: (nJ-1) x
     * nplist, row-major) */
    std::vector<realtype> s2llh;

    /** status code */
    int status = 0;

    /** total number of model parameters */
    const int np;
    /** number of fixed parameters */
    const int nk;
    /** number of states */
    const int nx;
    /** number of states with conservation laws applied */
    const int nx_solver;
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
    std::vector<ParameterScaling> pscale;
    /** flag indicating whether second order sensitivities were requested */
    const SecondOrderMode o2mode;
    /** sensitivity order */
    const SensitivityOrder sensi;
    /** sensitivity method */
    const SensitivityMethod sensi_meth;

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
