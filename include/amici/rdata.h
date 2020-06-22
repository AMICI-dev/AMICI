#ifndef AMICI_RDATA_H
#define AMICI_RDATA_H

#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/model.h"
#include "amici/misc.h"
#include "amici/forwardproblem.h"

#include <vector>

namespace amici {
class ReturnData;
class Solver;
class ExpData;
class ForwardProblem;
class BackwardProblem;
class SteadystateProblem;
} // namespace amici

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::ReturnData &u, unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief Stores all data to be returned by amici::runAmiciSimulation.
 *
 * NOTE: multidimensional arrays are stored in row-major order
 * (FORTRAN-style)
 */
class ReturnData {
  public:
    /**
     * @brief default constructor
     */
    ReturnData() = default;

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
     * @param nw see amici::Model::nw
     * @param pscale see amici::Model::pscale
     * @param o2mode see amici::Model::o2mode
     * @param sensi see amici::Solver::sensi
     * @param sensi_meth see amici::Solver::sensi_meth
     * @param rdrm see amici::Solver::rdata_reporting
     */
    ReturnData(std::vector<realtype> ts, int np, int nk, int nx, int nx_solver,
               int nxtrue, int ny, int nytrue, int nz, int nztrue, int ne,
               int nJ, int nplist, int nmaxevent, int nt, int newton_maxsteps,
               int nw,
               std::vector<ParameterScaling> pscale, SecondOrderMode o2mode,
               SensitivityOrder sensi, SensitivityMethod sensi_meth,
               RDataReporting rdrm);

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver instance
     * @param model model instance
     */
    ReturnData(Solver const &solver, const Model &model);

    ~ReturnData() = default;

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param preeq simulated preequilibration problem, pass nullptr to ignore
     * @param fwd simulated forward problem, pass nullptr to ignore
     * @param bwd simulat
     * ed backward problem, pass nullptr to ignore
     * @param posteq simulated postequilibration problem, pass nullptr to ignore
     * @param model matching model instance
     * @param solver matching solver instance
     * @param edata matching experimental data
     */
    void processSimulationObjects(SteadystateProblem const *preeq,
                                  ForwardProblem const *fwd,
                                  BackwardProblem const *bwd,
                                  SteadystateProblem const *posteq,
                                  Model &model, Solver const &solver,
                                  ExpData const *edata);

    /** timepoints (dimension: nt) */
    std::vector<realtype> ts;

    /** time derivative (dimension: nx) */
    std::vector<realtype> xdot;

    /**
     * Jacobian of differential equation right hand side (dimension: nx x nx,
     * row-major)
     */
    std::vector<realtype> J;

    /**
     * w data from the model (recurring terms in xdot, for imported SBML models
     * from python, this contains the flux vector)
     * (dimensions: nt x nw, row major)
     */
    std::vector<realtype> w;

    /** event output (dimension: nmaxevent x nz, row-major) */
    std::vector<realtype> z;

    /**
     * event output sigma standard deviation (dimension: nmaxevent x nz,
     * row-major)
     */
    std::vector<realtype> sigmaz;

    /**
     * parameter derivative of event output (dimension: nmaxevent x nz,
     * row-major)
     */
    std::vector<realtype> sz;

    /**
     * parameter derivative of event output standard deviation (dimension:
     * nmaxevent x nz, row-major)
     */
    std::vector<realtype> ssigmaz;

    /** event trigger output (dimension: nmaxevent x nz, row-major)*/
    std::vector<realtype> rz;

    /**
     * parameter derivative of event trigger output (dimension: nmaxevent x nz
     * x nplist, row-major)
     */
    std::vector<realtype> srz;

    /**
     * second order parameter derivative of event trigger output (dimension:
     * nmaxevent x nztrue x nplist x nplist, row-major)
     */
    std::vector<realtype> s2rz;

    /** state (dimension: nt x nx, row-major) */
    std::vector<realtype> x;

    /**
     * parameter derivative of state (dimension: nt x nplist x nx, row-major)
     */
    std::vector<realtype> sx;

    /** observable (dimension: nt x ny, row-major) */
    std::vector<realtype> y;

    /** observable standard deviation (dimension: nt x ny, row-major) */
    std::vector<realtype> sigmay;

    /**
     * parameter derivative of observable (dimension: nt x nplist x ny,
     * row-major)
     */
    std::vector<realtype> sy;

    /**
     * parameter derivative of observable standard deviation (dimension: nt x
     * nplist x ny, row-major)
     */
    std::vector<realtype> ssigmay;

    /** observable (dimension: nt*ny, row-major) */
    std::vector<realtype> res;

    /**
     * parameter derivative of residual (dimension: nt*ny x nplist, row-major)
     */
    std::vector<realtype> sres;

    /** fisher information matrix (dimension: nplist x nplist, row-major) */
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

    /**
     * number of linear solver convergence failures forward problem (dimension:
     * nt) */
    std::vector<int> numnonlinsolvconvfails;

    /**
     * number of linear solver convergence failures backwad problem (dimension:
     * nt) */
    std::vector<int> numnonlinsolvconvfailsB;

    /** employed order forward problem (dimension: nt) */
    std::vector<int> order;

    /** computation time of forward solve [ms] */
    double cpu_time = 0.0;

    /** computation time of backward solve [ms] */
    double cpu_timeB = 0.0;

    /** flags indicating success of steady state solver (preequilibration) */
    std::vector<SteadyStateStatus> preeq_status;

    /** computation time of the steady state solver [ms] (preequilibration) */
    double preeq_cpu_time = 0.0;

    /** computation time of the steady state solver of the backward problem [ms]
     *  (preequilibration) */
    double preeq_cpu_timeB = 0.0;

    /** flags indicating success of steady state solver  (postequilibration) */
    std::vector<SteadyStateStatus> posteq_status;

    /** computation time of the steady state solver [ms]  (postequilibration) */
    double posteq_cpu_time = 0.0;

    /** computation time of the steady state solverof the backward problem [ms]
     *  (postequilibration) */
    double posteq_cpu_timeB = 0.0;

    /**
     * number of Newton steps for steady state problem (preequilibration)
     * [newton, simulation, newton] (length = 3)
     */
    std::vector<int> preeq_numsteps;

    /**
     * number of linear steps by Newton step for steady state problem. this
     * will only be filled for iterative solvers (preequilibration)
     * (length = newton_maxsteps * 2)
     */
    std::vector<int> preeq_numlinsteps;

    /**
     * number of Newton steps for steady state problem (preequilibration)
     * [newton, simulation, newton] (length = 3) (postequilibration)
     */
    std::vector<int> posteq_numsteps;

    /**
     * number of linear steps by Newton step for steady state problem. this
     * will only be filled for iterative solvers (postequilibration)
     * (length = newton_maxsteps * 2)
     */
    std::vector<int> posteq_numlinsteps;

    /**
     * time at which steadystate was reached in the simulation based approach (preequilibration)
     */
    realtype preeq_t = NAN;

    /**
     * weighted root-mean-square of the rhs when steadystate
     * was reached (preequilibration)
     */
    realtype preeq_wrms = NAN;

    /**
     * time at which steadystate was reached in the simulation based approach (postequilibration)
     */
    realtype posteq_t = NAN;

    /**
     * weighted root-mean-square of the rhs when steadystate
     * was reached (postequilibration)
     */
    realtype posteq_wrms = NAN;

    /** initial state (dimension: nx) */
    std::vector<realtype> x0;

    /** preequilibration steady state found by Newton solver (dimension: nx) */
    std::vector<realtype> x_ss;

    /** initial sensitivities (dimension: nplist x nx, row-major) */
    std::vector<realtype> sx0;

    /**
     * preequilibration sensitivities found by Newton solver (dimension: nplist
     * x nx, row-major)
     */
    std::vector<realtype> sx_ss;

    /** loglikelihood value */
    realtype llh = 0.0;

    /** chi2 value */
    realtype chi2 = 0.0;

    /** parameter derivative of loglikelihood (dimension: nplist) */
    std::vector<realtype> sllh;

    /**
     * second order parameter derivative of loglikelihood (dimension: (nJ-1) x
     * nplist, row-major)
     */
    std::vector<realtype> s2llh;

    /** status code */
    int status = 0;

    /** total number of model parameters */
    int np{0};

    /** number of fixed parameters */
    int nk{0};

    /** number of states */
    int nx{0};

    /** number of states with conservation laws applied */
    int nx_solver{0};

    /** number of states in the unaugmented system */
    int nxtrue{0};

    /** number of observables */
    int ny{0};

    /** number of observables in the unaugmented system */
    int nytrue{0};

    /** number of event outputs */
    int nz{0};

    /** number of event outputs in the unaugmented system */
    int nztrue{0};

    /** number of events */
    int ne{0};

    /** dimension of the augmented objective function for 2nd order ASA */
    int nJ{0};

    /** number of parameter for which sensitivities were requested */
    int nplist{0};

    /** maximal number of occuring events (for every event type) */
    int nmaxevent{0};

    /** number of considered timepoints */
    int nt{0};

    /** number of columns in w */
    int nw{0};

    /** maximal number of newton iterations for steady state calculation */
    int newton_maxsteps{0};

    /** scaling of parameterization (lin,log,log10) */
    std::vector<ParameterScaling> pscale;

    /** flag indicating whether second order sensitivities were requested */
    SecondOrderMode o2mode{SecondOrderMode::none};

    /** sensitivity order */
    SensitivityOrder sensi{SensitivityOrder::none};

    /** sensitivity method */
    SensitivityMethod sensi_meth{SensitivityMethod::none};

    /** reporting mode */
    RDataReporting rdata_reporting{RDataReporting::full};

    /**
     * @brief Serialize ReturnData (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, ReturnData &r,
                                                unsigned int version);

  protected:

    /** timepoint for model evaluation*/
    realtype t;

    /** partial state vector, excluding states eliminated from conservation laws */
    AmiVector x_solver;

    /** partial time derivative of state vector, excluding states eliminated from conservation laws */
    AmiVector dx_solver;

    /** partial sensitivity state vector array, excluding states eliminated from
     * conservation laws */
    AmiVectorArray sx_solver;

    /** full state vector, including states eliminated from conservation laws */
    AmiVector x_rdata;

    /** full sensitivity state vector array, including states eliminated from
     * conservation laws */
    AmiVectorArray sx_rdata;

    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots;

    /**
     * @brief initializes storage for likelihood reporting mode
     */
    void initializeLikelihoodReporting();

    /**
     * @brief initializes storage for residual reporting mode
     */
    void initializeResidualReporting();

    /**
     * @brief initializes storage for full reporting mode
     */
    void initializeFullReporting();


    /**
     * @brief initialize values for chi2 and llh and derivatives
     */
    void initializeObjectiveFunction();

    /**
     * @brief extracts data from a preequilibration steadystateproblem
     * @param preeq Steadystateproblem for preequilibration
     * @param model Model instance to compute return values
     */
    void processPreEquilibration(SteadystateProblem const &preeq,
                                 Model &model);

    /**
     * @brief extracts data from a preequilibration steadystateproblem
     * @param posteq Steadystateproblem for postequilibration
     * @param model Model instance to compute return values
     * @param edata ExpData instance containing observable data
     */
    void processPostEquilibration(SteadystateProblem const &posteq,
                                  Model &model,
                                  ExpData const *edata);

    /**
     * @brief extracts results from forward problem
     * @param fwd forward problem
     * @param model model that was used for forward simulation
     * @param edata ExpData instance containing observable data
     */
    void processForwardProblem(ForwardProblem const &fwd,
                               Model &model,
                               ExpData const *edata);


    /**
     * @brief extracts results from backward problem
     * @param fwd forward problem
     * @param bwd backward problem
     * @param model model that was used for forward/backward simulation
     */
    void processBackwardProblem(ForwardProblem const &fwd,
                                BackwardProblem const &bwd,
                                Model &model);

    /**
     * @brief extracts results from solver
     * @param solver solver that was used for forward/backward simulation
     */
    void processSolver(Solver const &solver);

    /**
     * @brief Evaluates and stores the Jacobian and right hand side at final timepoint
     * @param problem forward problem or steadystate problem
     * @param model model that was used for forward/backward simulation
     */
    template <class T>
    void storeJacobianAndDerivativeInReturnData(T const &problem, Model &model)
    {
        readSimulationState(problem.getFinalSimulationState(), model);

        AmiVector xdot(nx_solver);
        if (!this->xdot.empty() || !this->J.empty())
            model.fxdot(t, x_solver, dx_solver, xdot);

        if (!this->xdot.empty())
            writeSlice(xdot, this->xdot);

        if (!this->J.empty()) {
            SUNMatrixWrapper J(nx_solver, nx_solver);
            model.fJ(t, 0.0, x_solver, dx_solver, xdot, J.get());
            // CVODES uses colmajor, so we need to transform to rowmajor
            for (int ix = 0; ix < model.nx_solver; ix++)
                for (int jx = 0; jx < model.nx_solver; jx++)
                    this->J.at(ix * model.nx_solver + jx) =
                        J.data()[ix + model.nx_solver * jx];
        }
    }
    /**
     * @brief sets member variables and model state according to provided simulation state
     * @param state simulation state provided by Problem
     * @param model model that was used for forward/backward simulation
     */
    void readSimulationState(SimulationState const &state, Model &model);

    /**
     * @brief Residual function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param edata ExpData instance containing observable data
     */
    void fres(int it, Model &model, const ExpData &edata);

    /**
     * @brief Chi-squared function
     * @param it time index
     */
    void fchi2(int it);

    /**
     * @brief Residual sensitivity function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param edata ExpData instance containing observable data
     */
    void fsres(int it, Model &model, const ExpData &edata);

    /**
     * @brief Fisher information matrix function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param edata ExpData instance containing observable data
     */
    void fFIM(int it, Model &model, const ExpData &edata);

    /**
     * @brief Set likelihood, state variables, outputs and respective
     * sensitivities to NaN (typically after integration failure)
     * @param it_start time index at which to start invalidating
     */
    void invalidate(int it_start);

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
     * @brief applies the chain rule to account for parameter transformation in
     * the sensitivities of simulation results
     * @param model Model from which the ReturnData was obtained
     */
    void applyChainRuleFactorToSimulationResults(const Model &model);


    /**
     * @brief Checks whether forward sensitivity analysis is performed
     * @return boolean indicator
     */
    bool computingFSA() const {
        return (sensi_meth == SensitivityMethod::forward &&
                sensi >= SensitivityOrder::first);
    }

    /**
     * @brief Extracts output information for data-points, expects that x_solver and sx_solver were
     * were set appropriately
     * @param it timepoint index
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getDataOutput(int it, Model &model, ExpData const *edata);

    /**
     * @brief Extracts data information for forward sensitivity analysis, expects that x_solver and
     * sx_solver were were set appropriately
     * @param it index of current timepoint
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getDataSensisFSA(int it, Model &model, ExpData const *edata);

    /**
     * @brief Extracts output information for events, expects that x_solver and sx_solver were
     * were set appropriately
     * @param iroot event index
     * @param t event timepoint
     * @param rootidx information about which roots fired (1 indicating fired, 0/-1 for not)
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getEventOutput(int iroot, realtype t, const std::vector<int> rootidx,
                        Model &model, ExpData const *edata);

    /**
     * @brief Extracts event information for forward sensitivity analysis, expects that x_solver and
     * sx_solver were set appropriately
     * @param iroot event index
     * @param ie index of event type
     * @param t event timepoint
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getEventSensisFSA(int iroot, int ie, realtype t, Model &model,
                           ExpData const *edata);
};

/**
 * @brief The ModelContext temporarily stores amici::Model::state
 * and restores it when going out of scope
 */
class ModelContext : public ContextManager {
  public:
    /**
     * @brief initialize backup of the original values.
     *
     * @param model
     */
    explicit ModelContext(Model *model);

    ModelContext &operator=(const ModelContext &other) = delete;

    ~ModelContext();

    /**
     * @brief Restore original state on constructor-supplied amici::Model.
     * Will be called during destruction. Explicit call is generally not
     * necessary.
     */
    void restore();

  private:
    Model *model = nullptr;
    ModelState original_state;
};


} // namespace amici

#endif /* _MY_RDATA */
