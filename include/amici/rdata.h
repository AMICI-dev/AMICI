#ifndef AMICI_RDATA_H
#define AMICI_RDATA_H

#include "amici/defines.h"
#include "amici/logging.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <vector>

namespace amici {
class ReturnData;
class Solver;
class ExpData;
class ForwardProblem;
class BackwardProblem;
class SteadystateProblem;
class SteadyStateBackwardProblem;
} // namespace amici

namespace boost::serialization {
template <class Archive>
void serialize(Archive& ar, amici::ReturnData& r, unsigned int version);
} // namespace boost::serialization

namespace amici {

/**
 * @brief Stores all data to be returned by amici::runAmiciSimulation.
 *
 * NOTE: multi-dimensional arrays are stored in row-major order (C-style)
 */
class ReturnData : public ModelDimensions {
  public:
    /**
     * @brief Default constructor
     */
    ReturnData() = default;

    /**
     * @brief Constructor
     * @param ts_ see amici::SimulationParameters::ts
     * @param model_dimensions_ Model dimensions
     * @param nmaxevent_ see amici::ModelDimensions::nmaxevent
     * @param newton_maxsteps_ see amici::Solver::newton_maxsteps
     * @param plist_ see amici::Model::getParameterList
     * @param pscale_ see amici::SimulationParameters::pscale
     * @param o2mode_ see amici::SimulationParameters::o2mode
     * @param sensi_ see amici::Solver::sensi
     * @param sensi_meth_ see amici::Solver::sensi_meth
     * @param rdrm_ see amici::Solver::rdata_reporting
     * @param quadratic_llh_ whether model defines a quadratic nllh and
     * computing res, sres and FIM makes sense
     * @param sigma_res_ indicates whether additional residuals are to be added
     * for each sigma
     * @param sigma_offset_ offset to ensure real-valuedness of sigma residuals
     */
    ReturnData(
        std::vector<realtype> ts_, ModelDimensions const& model_dimensions_,
        int nmaxevent_, int newton_maxsteps_, std::vector<int> plist_,
        std::vector<ParameterScaling> pscale_, SecondOrderMode o2mode_,
        SensitivityOrder sensi_, SensitivityMethod sensi_meth_,
        RDataReporting rdrm_, bool quadratic_llh_, bool sigma_res_,
        realtype sigma_offset_
    );

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver instance
     * @param model model instance
     */
    ReturnData(Solver const& solver, Model const& model);

    ~ReturnData() = default;

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param fwd simulated forward problem, pass `nullptr` to ignore
     * @param bwd simulated backward problem, pass `nullptr` to ignore
     * @param model matching model instance
     * @param solver matching solver instance
     * @param edata matching experimental data
     */
    void processSimulationObjects(
        ForwardProblem const* fwd, BackwardProblem const* bwd, Model& model,
        Solver const& solver, ExpData const* edata
    );
    /**
     * @brief Arbitrary (not necessarily unique) identifier.
     */
    std::string id;

    /**
     * Output or measurement timepoints (shape `nt`)
     */
    std::vector<realtype> ts;

    /** time derivative (shape `nx_solver`) evaluated at `t_last`. */
    std::vector<realtype> xdot;

    /**
     * @brief Jacobian of differential equation right hand side.
     *
     * The Jacobian of differential equation right hand side
     * (shape `nx_solver` x `nx_solver`, row-major) evaluated at `t_last`.
     *
     * The corresponding state variable IDs can be obtained from
     * `Model::getStateIdsSolver()`.
     */
    std::vector<realtype> J;

    /**
     * @brief Model expression values.
     *
     * Values of model expressions (recurring terms in xdot, for imported SBML
     * models from Python, this contains, e.g., the flux vector)
     * at timepoints `ReturnData::ts` (shape `nt` x `nw`, row major).
     *
     * The corresponding expression IDs can be obtained from
     * `Model::getExpressionIds()`.
     */
    std::vector<realtype> w;

    /** Event output (shape `nmaxevent` x `nz`, row-major) */
    std::vector<realtype> z;

    /**
     * Event output sigma standard deviation (shape `nmaxevent` x `nz`,
     * row-major)
     */
    std::vector<realtype> sigmaz;

    /**
     * Parameter derivative of event output
     * (shape `nmaxevent` x `nplist` x `nz`, row-major)
     */
    std::vector<realtype> sz;

    /**
     * Parameter derivative of event output standard deviation
     * (shape `nmaxevent` x `nplist` x `nz`, row-major)
     */
    std::vector<realtype> ssigmaz;

    /** Event trigger output (shape `nmaxevent` x `nz`, row-major)*/
    std::vector<realtype> rz;

    /**
     * Parameter derivative of event trigger output
     * (shape `nmaxevent` x `nplist` x `nz`, row-major)
     */
    std::vector<realtype> srz;

    /**
     * Second-order parameter derivative of event trigger output (shape
     * `nmaxevent` x `nztrue` x `nplist` x `nplist`, row-major)
     */
    std::vector<realtype> s2rz;

    /**
     * @brief Model state.
     *
     * The model state at timepoints `ReturnData::ts`
     * (shape `nt` x `nx_rdata`, row-major).
     *
     * The corresponding state variable IDs can be obtained from
     * `Model::getStateIds()`.
     */
    std::vector<realtype> x;

    /**
     * @brief State sensitivities.
     *
     * The derivative of the model state with respect to the chosen parameters
     * (see `Model::getParameterList()` or `ExpData::plist`)
     * at timepoints `ReturnData::ts`
     * (shape `nt` x `nplist` x `nx_rdata`, row-major).
     *
     * The corresponding state variable IDs can be obtained from
     * `Model::getStateIds()`.
     */
    std::vector<realtype> sx;

    /**
     * @brief Observables.
     *
     * The values of the observables at timepoints `ReturnData::ts`
     * (shape `nt` x `ny`, row-major).
     *
     * The corresponding observable IDs can be obtained from
     * `Model::getObservableIds()`.
     */
    std::vector<realtype> y;

    /** Observable standard deviation (shape `nt` x `ny`, row-major) */
    std::vector<realtype> sigmay;

    /**
     * @brief Observable sensitivities.
     *
     * The derivative of the observables with respect to the chosen parameters
     * (see `Model::getParameterList()` or `ExpData::plist`)
     * at timepoints `ReturnData::ts`
     * (shape `nt` x `nplist` x `ny`, row-major).
     *
     * The corresponding observable IDs can be obtained from
     * `Model::getObservableIds()`.
     */
    std::vector<realtype> sy;

    /**
     * Parameter derivative of observable standard deviation
     * (shape `nt` x `nplist` x `ny`, row-major)
     */
    std::vector<realtype> ssigmay;

    /** Residuals (shape `nt*ny`, row-major) */
    std::vector<realtype> res;

    /**
     * Parameter derivative of residual (shape `nt*ny` x `nplist`, row-major)
     */
    std::vector<realtype> sres;

    /** Fisher information matrix (shape `nplist` x `nplist`, row-major) */
    std::vector<realtype> FIM;

    /**
     * @brief Number of solver steps for the forward problem.
     *
     * Cumulative number of integration steps for the forward problem for each
     * output timepoint in `ReturnData::ts` (shape `nt`).
     */
    std::vector<int> numsteps;

    /**
     * @brief Number of solver steps for the backward problem.
     *
     * Cumulative number of integration steps for the backward problem for each
     * output timepoint in `ReturnData::ts` (shape `nt`).
     */
    std::vector<int> numstepsB;

    /** Number of right hand side evaluations forward problem (shape `nt`) */
    std::vector<int> numrhsevals;

    /** Number of right hand side evaluations backward problem (shape `nt`) */
    std::vector<int> numrhsevalsB;

    /** Number of error test failures forward problem (shape `nt`) */
    std::vector<int> numerrtestfails;

    /** Number of error test failures backward problem (shape `nt`) */
    std::vector<int> numerrtestfailsB;

    /**
     * Number of linear solver convergence failures forward problem (shape `nt`)
     */
    std::vector<int> numnonlinsolvconvfails;

    /**
     * Number of linear solver convergence failures backward problem (shape
     * `nt`)
     */
    std::vector<int> numnonlinsolvconvfailsB;

    /** Employed order forward problem (shape `nt`) */
    std::vector<int> order;

    /**
     * @brief Computation time of forward solve [ms]
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double cpu_time = 0.0;

    /**
     * @brief Computation time of backward solve [ms]
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double cpu_timeB = 0.0;

    /**
     * @brief Total CPU time from entering runAmiciSimulation until exiting [ms]
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double cpu_time_total = 0.0;

    /** Flags indicating success of steady-state solver (preequilibration) */
    std::vector<SteadyStateStatus> preeq_status;

    /**
     * @brief Computation time of the steady state solver [ms]
     * (pre-equilibration)
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double preeq_cpu_time = 0.0;

    /**
     * @brief Computation time of the steady state solver of the backward
     * problem [ms] (pre-equilibration)
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double preeq_cpu_timeB = 0.0;

    /** Flags indicating success of steady-state solver  (postequilibration) */
    std::vector<SteadyStateStatus> posteq_status;

    /**
     * @brief Computation time of the steady-state solver [ms]
     * (postequilibration)
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double posteq_cpu_time = 0.0;

    /**
     * @brief Computation time of the steady-state solver of the backward
     * problem [ms] (postequilibration)
     *
     * .. warning::
     *      If AMICI was built without boost, this tracks the CPU-time of the
     *      current process. Therefore, in a multi-threaded context, this value
     *      may be incorrect.
     *
     */
    double posteq_cpu_timeB = 0.0;

    /**
     * @brief Number of Newton steps for pre-equilibration.
     *
     * [newton, simulation, newton] (length = 3)
     */
    std::vector<int> preeq_numsteps;

    /**
     * Number of simulation steps for adjoint pre-equilibration problem
     * [== 0 if analytical solution worked, > 0 otherwise]
     */
    int preeq_numstepsB = 0;

    /**
     * Number of Newton steps for post-equilibration
     * [newton, simulation, newton].
     */
    std::vector<int> posteq_numsteps;

    /**
     * Number of simulation steps for the post-equilibration backward simulation
     * [== 0 if analytical solution worked, > 0 otherwise]
     */
    int posteq_numstepsB = 0;

    /**
     * Model time at which the pre-equilibration steady state was reached via
     * simulation.
     */
    realtype preeq_t = NAN;

    /**
     * Weighted root-mean-square of the rhs at pre-equilibration steady state.
     */
    realtype preeq_wrms = NAN;

    /**
     * Model time at which the post-equilibration steady state was reached via
     * simulation.
     */
    realtype posteq_t = NAN;

    /**
     * Weighted root-mean-square of the rhs at post-equilibration steady state.
     */
    realtype posteq_wrms = NAN;

    /**
     * @brief Initial state of the main simulation (shape `nx_rdata`).
     *
     * The corresponding state variable IDs can be obtained from
     * `Model::getStateIds()`.
     */
    std::vector<realtype> x0;

    /**
     * @brief Pre-equilibration steady state.
     *
     * The values of the state variables at the pre-equilibration steady state
     * (shape `nx_rdata`).
     * The corresponding state variable IDs can be obtained from
     * `Model::getStateIds()`.
     */
    std::vector<realtype> x_ss;

    /**
     * @brief Initial state sensitivities for the main simulation.
     *
     * (shape `nplist` x `nx_rdata`, row-major).
     */
    std::vector<realtype> sx0;

    /**
     * @brief Pre-equilibration steady state sensitivities.

     * Sensitivities of the pre-equilibration steady state with respect to the
     * selected parameters.
     * (shape `nplist` x `nx_rdata`, row-major)
     */
    std::vector<realtype> sx_ss;

    /** Log-likelihood value */
    realtype llh = 0.0;

    /** \f$\chi^2\f$ value */
    realtype chi2 = 0.0;

    /** Parameter derivative of log-likelihood (shape `nplist`) */
    std::vector<realtype> sllh;

    /**
     * Second-order parameter derivative of log-likelihood
     * (shape `nJ-1` x `nplist`, row-major)
     */
    std::vector<realtype> s2llh;

    /**
     * @brief Simulation status code.
     *
     * One of:
     *
     * * AMICI_SUCCESS, indicating successful simulation
     * * AMICI_MAX_TIME_EXCEEDED, indicating that the simulation did not finish
     *   within the allowed time (see Solver.{set,get}MaxTime)
     * * AMICI_ERROR, indicating that some error occurred during simulation
     *   (a more detailed error message will have been printed).
     * * AMICI_NOT_RUN, if no simulation was started
     */
    int status = AMICI_NOT_RUN;

    /**
     * @brief Number of state variables.
     *
     * (alias `nx_rdata`, kept for backward compatibility)
     */
    int nx{0};

    /**
     * Number of state variables in the unaugmented system
     * (alias nxtrue_rdata, kept for backward compatibility)
     */
    int nxtrue{0};

    /** Number of parameters w.r.t. which sensitivities were requested */
    int nplist{0};

    /** Maximal number of occurring events (for every event type) */
    int nmaxevent{0};

    /** Number of output timepoints (length of `ReturnData::ts`). */
    int nt{0};

    /** maximal number of newton iterations for steady state calculation */
    int newton_maxsteps{0};

    /** Scaling of model parameters. */
    std::vector<ParameterScaling> pscale;

    /** Flag indicating whether second-order sensitivities were requested. */
    SecondOrderMode o2mode{SecondOrderMode::none};

    /** Sensitivity order. */
    SensitivityOrder sensi{SensitivityOrder::none};

    /** Sensitivity method. */
    SensitivityMethod sensi_meth{SensitivityMethod::none};

    /** Reporting mode. */
    RDataReporting rdata_reporting{RDataReporting::full};

    /**
     * @brief Serialize ReturnData (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(
        Archive& ar, ReturnData& r, unsigned int version
    );

    /**
     * Boolean indicating whether residuals for standard deviations have been
     * added.
     */
    bool sigma_res{false};

    /** Log messages. */
    std::vector<LogItem> messages;

    /** The final internal time of the solver. */
    realtype t_last{std::numeric_limits<realtype>::quiet_NaN()};

    /**
     * @brief Indices of the parameters w.r.t. which sensitivities were
     * computed.
     *
     * The indices refer to parameter IDs in Model::getParameterIds().
     */
    std::vector<int> plist;

  protected:
    /** offset for sigma_residuals */
    realtype sigma_offset{0.0};

    /** array of number of found roots for a certain event type
     * (shape `ne`) */
    std::vector<int> nroots_;

    /**
     * @brief initializes storage for likelihood reporting mode
     * @param quadratic_llh whether model defines a quadratic nllh and computing
     * res, sres and FIM makes sense.
     */
    void initializeLikelihoodReporting(bool quadratic_llh);

    /**
     * @brief initializes storage for observables + likelihood reporting mode
     * @param quadratic_llh whether model defines a quadratic nllh and computing
     * res, sres and FIM makes sense.
     */
    void initializeObservablesLikelihoodReporting(bool quadratic_llh);

    /**
     * @brief initializes storage for residual reporting mode
     * @param enable_res whether residuals are to be computed
     */
    void initializeResidualReporting(bool enable_res);

    /**
     * @brief initializes storage for full reporting mode
     * @param enable_fim whether FIM Hessian approximation is to be computed
     */
    void initializeFullReporting(bool enable_fim);

    /**
     * @brief initialize values for chi2 and llh and derivatives
     * @param enable_chi2 whether chi2 values are to be computed
     */
    void initializeObjectiveFunction(bool enable_chi2);

    /**
     * @brief extracts data from a preequilibration SteadystateProblem
     * @param preeq SteadystateProblem for preequilibration
     * @param preeq_bwd SteadyStateBackwardProblem from preequilibration
     * @param model Model instance to compute return values
     */
    void processPreEquilibration(
        SteadystateProblem const& preeq,
        SteadyStateBackwardProblem const* preeq_bwd, Model& model
    );

    /**
     * @brief extracts data from a preequilibration SteadystateProblem
     * @param posteq SteadystateProblem for postequilibration
     * @param posteq_bwd SteadyStateBackwardProblem from postequilibration
     * @param model Model instance to compute return values
     * @param edata ExpData instance containing observable data
     */
    void processPostEquilibration(
        SteadystateProblem const& posteq,
        SteadyStateBackwardProblem const* posteq_bwd, Model& model,
        ExpData const* edata
    );

    /**
     * @brief extracts results from forward problem
     * @param fwd forward problem
     * @param model model that was used for forward simulation
     * @param edata ExpData instance containing observable data
     */
    void processForwardProblem(
        ForwardProblem const& fwd, Model& model, ExpData const* edata
    );

    /**
     * @brief extracts results from backward problem
     * @param fwd forward problem
     * @param bwd backward problem
     * @param preeq SteadyStateProblem for preequilibration
     * @param preeq_bwd SteadyStateBackwardProblem for preequilibration
     * @param model model that was used for forward/backward simulation
     */
    void processBackwardProblem(
        ForwardProblem const& fwd, BackwardProblem const& bwd,
        SteadystateProblem const* preeq,
        SteadyStateBackwardProblem const* preeq_bwd, Model& model
    );

    /**
     * @brief extracts results from solver
     * @param solver solver that was used for forward/backward simulation
     */
    void processSolver(Solver const& solver);

    /**
     * @brief Evaluates and stores the Jacobian and right hand side at final
     * timepoint
     * @param problem forward problem or steadystate problem
     * @param model model that was used for forward/backward simulation
     */
    template <class T>
    void
    storeJacobianAndDerivativeInReturnData(T const& problem, Model& model) {
        auto const& simulation_state = problem.getFinalSimulationState();
        model.setModelState(simulation_state.mod);
        auto const& sol = simulation_state.sol;
        sundials::Context sunctx;
        AmiVector xdot(nx_solver, sunctx);
        if (!this->xdot.empty() || !this->J.empty())
            model.fxdot(sol.t, sol.x, sol.dx, xdot);

        if (!this->xdot.empty())
            writeSlice(xdot, this->xdot);

        if (!this->J.empty()) {
            SUNMatrixWrapper J(nx_solver, nx_solver, sunctx);
            model.fJ(sol.t, 0.0, sol.x, sol.dx, xdot, J);
            // CVODES uses colmajor, so we need to transform to rowmajor
            for (int ix = 0; ix < model.nx_solver; ix++)
                for (int jx = 0; jx < model.nx_solver; jx++)
                    this->J.at(ix * model.nx_solver + jx)
                        = J.data()[ix + model.nx_solver * jx];
        }
    }

    /**
     * @brief Residual function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param sol Solution state the timepoint `it`
     * @param edata ExpData instance containing observable data
     */
    void
    fres(int it, Model& model, SolutionState const& sol, ExpData const& edata);

    /**
     * @brief Chi-squared function
     * @param it time index
     * @param edata ExpData instance containing observable data
     */
    void fchi2(int it, ExpData const& edata);

    /**
     * @brief Residual sensitivity function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param sol solution state the timepoint `it`
     * @param edata ExpData instance containing observable data
     */
    void
    fsres(int it, Model& model, SolutionState const& sol, ExpData const& edata);

    /**
     * @brief Fisher information matrix function
     * @param it time index
     * @param model model that was used for forward/backward simulation
     * @param sol Solution state the timepoint `it`
     * @param edata ExpData instance containing observable data
     */
    void
    fFIM(int it, Model& model, SolutionState const& sol, ExpData const& edata);

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
    void applyChainRuleFactorToSimulationResults(Model const& model);

    /**
     * @brief Checks whether forward sensitivity analysis is performed
     * @return boolean indicator
     */
    [[nodiscard]] bool computingFSA() const {
        return (
            sensi_meth == SensitivityMethod::forward
            && sensi >= SensitivityOrder::first
        );
    }

    /**
     * @brief Extracts output information for data-points, expects that
     * the model state was set appropriately
     * @param it timepoint index
     * @param model model that was used in forward solve
     * @param sol solution state the timepoint `it`
     * @param edata ExpData instance carrying experimental data
     */
    void getDataOutput(
        int it, Model& model, SolutionState const& sol, ExpData const* edata
    );

    /**
     * @brief Extracts data information for forward sensitivity analysis,
     * expects that the model state was set appropriately
     * @param it index of current timepoint
     * @param model model that was used in forward solve
     * @param sol Solution state the timepoint `it`
     * @param edata ExpData instance carrying experimental data
     */
    void getDataSensisFSA(
        int it, Model& model, SolutionState const& sol, ExpData const* edata
    );

    /**
     * @brief Extracts output information for events, expects that the model
     * state was set appropriately
     * @param rootidx information about which roots fired
     * (1 indicating fired, 0/-1 for not)
     * @param model model that was used in forward solve
     * @param sol Solution state the timepoint `it`
     * @param edata ExpData instance carrying experimental data
     */
    void getEventOutput(
        std::vector<int> const& rootidx, Model& model, SolutionState const& sol,
        ExpData const* edata
    );

    /**
     * @brief Extracts event information for forward sensitivity analysis,
     * expects the model state was set appropriately
     * @param ie index of event type
     * @param model model that was used in forward solve
     * @param sol Solution state the timepoint `it`
     * @param edata ExpData instance carrying experimental data
     */
    void getEventSensisFSA(
        int ie, Model& model, SolutionState const& sol, ExpData const* edata
    );

    /**
     * @brief Updates contribution to likelihood from quadratures (xQB),
     * if preequilibration was run in adjoint mode
     * @param model model that was used for forward/backward simulation
     * @param sx0 State sensitivities at pre-equilibration steady state
     * @param xB Adjoint state from pre-equilibration.
     * @param llhS0 contribution to likelihood for initial state sensitivities
     * of preequilibration
     */
    void handleSx0Backward(
        Model const& model, AmiVectorArray const& sx0, AmiVector const& xB,
        std::vector<realtype>& llhS0
    ) const;

    /**
     * @brief Updates contribution to likelihood for initial state sensitivities
     * (llhS0), if no preequilibration was run or if forward sensitivities were
     * used
     * @param model model that was used for forward/backward simulation
     * @param sol Solution state the timepoint `it`
     * @param llhS0 contribution to likelihood for initial state sensitivities
     * @param xB vector with final adjoint state
     * (excluding conservation laws)
     */
    void handleSx0Forward(
        Model const& model, SolutionState const& sol,
        std::vector<realtype>& llhS0, AmiVector const& xB
    ) const;
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
    explicit ModelContext(Model* model);

    ModelContext& operator=(ModelContext const& other) = delete;

    ~ModelContext() noexcept(false);

    /**
     * @brief Restore original state on constructor-supplied amici::Model.
     * Will be called during destruction. Explicit call is generally not
     * necessary.
     */
    void restore();

  private:
    Model* model_{nullptr};
    ModelState original_state_;
};

} // namespace amici

#endif // AMICI_RDATA_H
