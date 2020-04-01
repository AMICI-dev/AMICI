#ifndef AMICI_RDATA_H
#define AMICI_RDATA_H

#include "amici/defines.h"
#include "amici/vector.h"


#include <vector>

namespace amici {
class Model;
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
     */
    ReturnData(std::vector<realtype> ts, int np, int nk, int nx, int nx_solver,
               int nxtrue, int ny, int nytrue, int nz, int nztrue, int ne,
               int nJ, int nplist, int nmaxevent, int nt, int newton_maxsteps,
               int nw,
               std::vector<ParameterScaling> pscale, SecondOrderMode o2mode,
               SensitivityOrder sensi, SensitivityMethod sensi_meth);

    /**
     * @brief constructor that uses information from model and solver to
     * appropriately initialize fields
     * @param solver solver instance
     * @param model model instance
     */
    ReturnData(Solver const &solver, const Model &model);

    ~ReturnData() = default;

    /**
     * @brief initializeObjectiveFunction
     */
    void initializeObjectiveFunction();

    /**
     * @brief extracts data from a preequilibration steadystateproblem
     * @param preeq Steadystateproblem
     * @param model Model from which the ReturnData was obtained
     */
    void processPreequilibration(SteadystateProblem const *preeq, Model *model);

    /**
     * @brief extracts results from forward problem
     * @param fwd forward problem
     * @param model model that was used for forward simulation
     * @param edata ExpData instance containing observable data
     */
    void processForwardProblem(ForwardProblem const &fwd,
                               Model *model,
                               ExpData const *edata);
    
    
    /**
     * @brief extracts results from backward problem
     * @param fwd forward problem
     * @param bwd backward problem
     * @param model model that was used for forward/backward simulation
     */
    void processBackwardProblem(ForwardProblem const &fwd,
                                BackwardProblem const &bwd,
                                Model *model);
    
    /**
     * @brief extracts results from solver
     * @param solver solverl that was used for forward/backward simulation
     */
    void processSolver(Solver const &solver);
    
    /**
     * @brief Evaluates the Jacobian and differential equation right hand side,
     * stores it in RetunData
     * @param solver solverl that was used for forward/backward simulation
     * @param model model that was used for forward/backward simulation
     */
    void storeJacobianAndDerivativeInReturnData(Solver const &solver,
                                                Model *model);
    
    
    /**
     * @brief applies the chain rule to account for parameter transformation in
     * the sensitivities of simulation results
     * @param model Model from which the ReturnData was obtained
     */
    void applyChainRuleFactorToSimulationResults(const Model *model);

    /**
     * @brief Residual function
     * @param it time index
     * @param edata ExpData instance containing observable data
     */
    void fres(int it, const ExpData &edata);

    /**
     * @brief Chi-squared function
     * @param it time index
     */
    void fchi2(int it);

    /**
     * @brief Residual sensitivity function
     * @param it time index
     * @param edata ExpData instance containing observable data
     */
    void fsres(int it, const ExpData &edata);

    /**
     * @brief Fisher information matrix function
     * @param it time index
     */
    void fFIM(int it);

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

    /** flag indicating success of Newton solver */
    int newton_status = 0;

    /** computation time of the Newton solver [ms] */
    double newton_cpu_time = 0.0;

    /**
     * number of Newton steps for steady state problem
     * [newton, simulation, newton] (length = 3)
     */
    std::vector<int> newton_numsteps;

    /**
     * number of linear steps by Newton step for steady state problem. this
     * will only be filled for iterative solvers (length = newton_maxsteps * 2)
     */
    std::vector<int> newton_numlinsteps;

    /**
     * time at which steadystate was reached in the simulation based approach
     */
    realtype t_steadystate = NAN;

    /**
     * weighted root-mean-square of the rhs when steadystate
     * was reached
     */
    realtype wrms_steadystate = NAN;

    /**
     * weighted root-mean-square of the rhs when steadystate
     * was reached
     */
    realtype wrms_sensi_steadystate = NAN;

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
    /** state vector, including states eliminated from conservation laws
     * (dimension: nx) */
    AmiVector x_rdata;
    
    /** full sensitivity state vector array, including states eliminated from
     * conservation laws (dimension: nx x nplist, row-major) */
    AmiVectorArray sx_rdata;
    
    /** array of number of found roots for a certain event type
     * (dimension: ne) */
    std::vector<int> nroots;
    


    /**
     * @brief Extracts data information for forward sensitivity analysis
     * @param it index of current timepoint
     * @param fwd forward problem
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getDataSensisFSA(int it, ForwardProblem const &fwd, Model *model,
                          ExpData const *edata);
    
    /**
     * @brief Extracts output information for data-points
     * @param it timepoint index
     * @param fwd forward problem
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getDataOutput(int it, ForwardProblem const &fwd, Model *model,
                       ExpData const *edata);
    
    /**
     * @brief Extracts output information for events
     * @param iroot event index
     * @param fwd forward problem
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getEventOutput(int iroot, ForwardProblem const &fwd, Model *model,
                        ExpData const *edata);
    
    /**
     * @brief Extracts event information for forward sensitivity analysis
     *
     * @param iroot event index
     * @param ie index of event type
     * @param fwd forward problem
     * @param model model that was used in forward solve
     * @param edata ExpData instance carrying experimental data
     */
    void getEventSensisFSA(int iroot, int ie, ForwardProblem const &fwd,
                           Model *model, ExpData const *edata);
};

} // namespace amici

#endif /* _MY_RDATA */
