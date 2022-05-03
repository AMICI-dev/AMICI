#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"

#include "sundials/sundials_types.h"

namespace amici {

class Solver;

namespace model_TPL_MODELNAME {

extern std::array<const char*, TPL_NP> parameterNames;
extern std::array<const char*, TPL_NK> fixedParameterNames;
extern std::array<const char*, TPL_NX_RDATA> stateNames;
extern std::array<const char*, TPL_NY> observableNames;
extern std::array<const ObservableScaling, TPL_NY> observableScalings;
extern std::array<const char*, TPL_NW> expressionNames;
extern std::array<const char*, TPL_NP> parameterIds;
extern std::array<const char*, TPL_NK> fixedParameterIds;
extern std::array<const char*, TPL_NX_RDATA> stateIds;
extern std::array<const char*, TPL_NY> observableIds;
extern std::array<const char*, TPL_NW> expressionIds;
extern std::array<int, TPL_NX_SOLVER> stateIdxsSolver;
extern std::array<bool, TPL_NEVENT> rootInitialValues;

TPL_JY_DEF
TPL_DJYDSIGMA_DEF
TPL_DJYDY_DEF
TPL_DJYDY_COLPTRS_DEF
TPL_DJYDY_ROWVALS_DEF
TPL_ROOT_DEF
TPL_DWDP_DEF
TPL_DWDP_COLPTRS_DEF
TPL_DWDP_ROWVALS_DEF
TPL_DWDX_DEF
TPL_DWDX_COLPTRS_DEF
TPL_DWDX_ROWVALS_DEF
TPL_DWDW_DEF
TPL_DWDW_COLPTRS_DEF
TPL_DWDW_ROWVALS_DEF
TPL_DXDOTDW_DEF
TPL_DXDOTDW_COLPTRS_DEF
TPL_DXDOTDW_ROWVALS_DEF
TPL_DXDOTDP_EXPLICIT_DEF
TPL_DXDOTDP_EXPLICIT_COLPTRS_DEF
TPL_DXDOTDP_EXPLICIT_ROWVALS_DEF
TPL_DXDOTDX_EXPLICIT_DEF
TPL_DXDOTDX_EXPLICIT_COLPTRS_DEF
TPL_DXDOTDX_EXPLICIT_ROWVALS_DEF
TPL_DYDX_DEF
TPL_DYDP_DEF
TPL_SIGMAY_DEF
TPL_DSIGMAYDP_DEF
TPL_DSIGMAYDY_DEF
TPL_W_DEF
TPL_X0_DEF
TPL_X0_FIXEDPARAMETERS_DEF
TPL_SX0_DEF
TPL_SX0_FIXEDPARAMETERS_DEF
TPL_XDOT_DEF
TPL_Y_DEF
TPL_STAU_DEF
TPL_DELTAX_DEF
TPL_DELTASX_DEF
TPL_X_RDATA_DEF
TPL_X_SOLVER_DEF
TPL_TOTAL_CL_DEF
TPL_DX_RDATADX_SOLVER_DEF
TPL_DX_RDATADX_SOLVER_COLPTRS_DEF
TPL_DX_RDATADX_SOLVER_ROWVALS_DEF
TPL_DX_RDATADP_DEF
TPL_DX_RDATADTCL_DEF
TPL_DX_RDATADTCL_COLPTRS_DEF
TPL_DX_RDATADTCL_ROWVALS_DEF
TPL_DTOTAL_CLDP_DEF
TPL_DTOTAL_CLDX_RDATA_DEF
TPL_DTOTAL_CLDX_RDATA_COLPTRS_DEF
TPL_DTOTAL_CLDX_RDATA_ROWVALS_DEF
/**
 * @brief AMICI-generated model subclass.
 */
class Model_TPL_MODELNAME : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_TPL_MODELNAME()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  TPL_NX_RDATA,                            // nx_rdata
                  TPL_NXTRUE_RDATA,                        // nxtrue_rdata
                  TPL_NX_SOLVER,                           // nx_solver
                  TPL_NXTRUE_SOLVER,                       // nxtrue_solver
                  TPL_NX_SOLVER_REINIT,                    // nx_solver_reinit
                  TPL_NP,                                  // np
                  TPL_NK,                                  // nk
                  TPL_NY,                                  // ny
                  TPL_NYTRUE,                              // nytrue
                  TPL_NZ,                                  // nz
                  TPL_NZTRUE,                              // nztrue
                  TPL_NEVENT,                              // nevent
                  TPL_NOBJECTIVE,                          // nobjective
                  TPL_NW,                                  // nw
                  TPL_NDWDX,                               // ndwdx
                  TPL_NDWDP,                               // ndwdp
                  TPL_NDWDW,                               // ndwdw
                  TPL_NDXDOTDW,                            // ndxdotdw
                  TPL_NDJYDY,                              // ndjydy
                  TPL_NDXRDATADXSOLVER,                    // ndxrdatadxsolver
                  TPL_NDXRDATADTCL,                        // ndxrdatadtcl
                  TPL_NDTOTALCLDXRDATA,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  TPL_UBW,                                 // ubw
                  TPL_LBW                                  // lbw
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{TPL_FIXED_PARAMETERS}, // fixedParameters
                  std::vector<realtype>{TPL_PARAMETERS}        // dynamic parameters
              ),
              TPL_O2MODE,                                  // o2mode
              std::vector<realtype>(TPL_NX_SOLVER, 0.0),   // idlist
              std::vector<int>{},                          // z2event
              true,                                        // pythonGenerated
              TPL_NDXDOTDP_EXPLICIT,                       // ndxdotdp_explicit
              TPL_NDXDOTDX_EXPLICIT,                       // ndxdotdx_explicit
              TPL_W_RECURSION_DEPTH                        // w_recursion_depth
          ) {
                 root_initial_values_ = std::vector<bool>(
                     rootInitialValues.begin(), rootInitialValues.end()
                 );
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_TPL_MODELNAME(*this);
    }

    /**
     * @brief model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fJrz(realtype *nllh, const int iz, const realtype *p,
              const realtype *k, const realtype *rz,
              const realtype *sigmaz) override {}

    TPL_JY_IMPL

    /**
     * @brief model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     */
    void fJz(realtype *nllh, const int iz, const realtype *p,
             const realtype *k, const realtype *z,
             const realtype *sigmaz, const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fdJrzdsigma(realtype *dJrzdsigma, const int iz,
                     const realtype *p, const realtype *k,
                     const realtype *rz,
                     const realtype *sigmaz) override {}

    /**
     * @brief model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     */
    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p,
                 const realtype *k, const realtype *rz,
                 const realtype *sigmaz) override {}

    TPL_DJYDSIGMA_IMPL

    /**
     * @brief model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     */
    void fdJzdsigma(realtype *dJzdsigma, const int iz,
                    const realtype *p, const realtype *k,
                    const realtype *z, const realtype *sigmaz,
                    const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative
     * log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     */
    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p,
                const realtype *k, const realtype *z,
                const realtype *sigmaz, const realtype *mz) override {}

    /**
     * @brief model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     */
    void fdeltaqB(realtype *deltaqB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ip,
                  const int ie, const realtype *xdot,
                  const realtype *xdot_old,
                  const realtype *xB) override {}

    TPL_DELTASX_IMPL

    TPL_DELTAX_IMPL

    /**
     * @brief model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     */
    void fdeltaxB(realtype *deltaxB, const realtype t,
                  const realtype *x, const realtype *p,
                  const realtype *k, const realtype *h, const int ie,
                  const realtype *xdot, const realtype *xdot_old,
                  const realtype *xB) override {}

    /**
     * @brief model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters
     * p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     */
    void fdrzdp(realtype *drzdp, const int ie, const realtype t,
                const realtype *x, const realtype *p, const realtype *k,
                const realtype *h, const int ip) override {}

    /**
     * @brief model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fdrzdx(realtype *drzdx, const int ie, const realtype t,
                const realtype *x, const realtype *p, const realtype *k,
                const realtype *h) override {}

    TPL_DSIGMAYDP_IMPL

    TPL_DSIGMAYDY_IMPL

    /**
     * @brief model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event
     * measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     */
    void fdsigmazdp(realtype *dsigmazdp, const realtype t,
                    const realtype *p, const realtype *k,
                    const int ip) override {}

    TPL_DJYDY_IMPL
    TPL_DJYDY_COLPTRS_IMPL
    TPL_DJYDY_ROWVALS_IMPL

    TPL_DWDP_IMPL
    TPL_DWDP_COLPTRS_IMPL
    TPL_DWDP_ROWVALS_IMPL

    TPL_DWDX_IMPL
    TPL_DWDX_COLPTRS_IMPL
    TPL_DWDX_ROWVALS_IMPL

    TPL_DWDW_IMPL
    TPL_DWDW_COLPTRS_IMPL
    TPL_DWDW_ROWVALS_IMPL

    TPL_DXDOTDW_IMPL
    TPL_DXDOTDW_COLPTRS_IMPL
    TPL_DXDOTDW_ROWVALS_IMPL

    TPL_DXDOTDP_EXPLICIT_IMPL
    TPL_DXDOTDP_EXPLICIT_COLPTRS_IMPL
    TPL_DXDOTDP_EXPLICIT_ROWVALS_IMPL

    TPL_DXDOTDX_EXPLICIT_IMPL
    TPL_DXDOTDX_EXPLICIT_COLPTRS_IMPL
    TPL_DXDOTDX_EXPLICIT_ROWVALS_IMPL

    TPL_DYDX_IMPL

    TPL_DYDP_IMPL

    /**
     * @brief model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model
     * parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     */
    void fdzdp(realtype *dzdp, const int ie, const realtype t,
               const realtype *x, const realtype *p, const realtype *k,
               const realtype *h, const int ip) override {}

    /**
     * @brief model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model
     * states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fdzdx(realtype *dzdx, const int ie, const realtype t,
               const realtype *x, const realtype *p, const realtype *k,
               const realtype *h) override {}

    TPL_ROOT_IMPL

    /**
     * @brief model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events
     * not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void frz(realtype *rz, const int ie, const realtype t,
             const realtype *x, const realtype *p, const realtype *k,
             const realtype *h) override {}

    TPL_SIGMAY_IMPL

    /**
     * @brief model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     */
    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p,
                 const realtype *k) override {}

    /**
     * @brief model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heaviside vector
     * @param ip sensitivity index
     */
    void fsrz(realtype *srz, const int ie, const realtype t,
              const realtype *x, const realtype *p, const realtype *k,
              const realtype *h, const realtype *sx,
              const int ip) override {}

    TPL_STAU_IMPL
    TPL_SX0_IMPL
    TPL_SX0_FIXEDPARAMETERS_IMPL

    /**
     * @brief model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     */
    void fsz(realtype *sz, const int ie, const realtype t,
             const realtype *x, const realtype *p, const realtype *k,
             const realtype *h, const realtype *sx,
             const int ip) override {}

    TPL_W_IMPL

    TPL_X0_IMPL

    TPL_X0_FIXEDPARAMETERS_IMPL

    TPL_XDOT_IMPL

    TPL_Y_IMPL

    /**
     * @brief model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heaviside vector
     */
    void fz(realtype *z, const int ie, const realtype t,
            const realtype *x, const realtype *p, const realtype *k,
            const realtype *h) override {}

    TPL_X_RDATA_IMPL

    TPL_X_SOLVER_IMPL

    TPL_TOTAL_CL_IMPL

    TPL_DX_RDATADX_SOLVER_IMPL
    TPL_DX_RDATADX_SOLVER_COLPTRS_IMPL
    TPL_DX_RDATADX_SOLVER_ROWVALS_IMPL

    TPL_DX_RDATADP_IMPL

    TPL_DX_RDATADTCL_IMPL
    TPL_DX_RDATADTCL_COLPTRS_IMPL
    TPL_DX_RDATADTCL_ROWVALS_IMPL

    TPL_DTOTAL_CLDP_IMPL

    TPL_DTOTAL_CLDX_RDATA_IMPL
    TPL_DTOTAL_CLDX_RDATA_COLPTRS_IMPL
    TPL_DTOTAL_CLDX_RDATA_ROWVALS_IMPL

    std::string getName() const override {
        return "TPL_MODELNAME";
    }

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    std::vector<std::string> getParameterNames() const override {
        return std::vector<std::string>(parameterNames.begin(),
                                        parameterNames.end());
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    std::vector<std::string> getStateNames() const override {
        return std::vector<std::string>(stateNames.begin(), stateNames.end());
    }

    /**
     * @brief Get names of the solver states
     * @return the names
     */
    std::vector<std::string> getStateNamesSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateNames[idx]);
        }
        return result;
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    std::vector<std::string> getFixedParameterNames() const override {
        return std::vector<std::string>(fixedParameterNames.begin(),
                                        fixedParameterNames.end());
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    std::vector<std::string> getObservableNames() const override {
        return std::vector<std::string>(observableNames.begin(),
                                        observableNames.end());
    }

    /**
     * @brief Get names of model expressions
     * @return Expression names
     */
    std::vector<std::string> getExpressionNames() const override {
        return std::vector<std::string>(expressionNames.begin(),
                                        expressionNames.end());
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    std::vector<std::string> getParameterIds() const override {
        return std::vector<std::string>(parameterIds.begin(),
                                        parameterIds.end());
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    std::vector<std::string> getStateIds() const override {
        return std::vector<std::string>(stateIds.begin(), stateIds.end());
    }

    /**
     * @brief Get ids of the solver states
     * @return the ids
     */
    std::vector<std::string> getStateIdsSolver() const override {
        std::vector<std::string> result;
        result.reserve(stateIdxsSolver.size());
        for(auto idx: stateIdxsSolver) {
            result.push_back(stateIds[idx]);
        }
        return result;
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    std::vector<std::string> getFixedParameterIds() const override {
        return std::vector<std::string>(fixedParameterIds.begin(),
                                        fixedParameterIds.end());
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    std::vector<std::string> getObservableIds() const override {
        return std::vector<std::string>(observableIds.begin(),
                                        observableIds.end());
    }

    /**
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    std::vector<std::string> getExpressionIds() const override {
        return std::vector<std::string>(expressionIds.begin(),
                                        expressionIds.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending
     * on fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    bool isFixedParameterStateReinitializationAllowed() const override {
        return TPL_REINIT_FIXPAR_INITCOND;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string getAmiciVersion() const override {
        return "TPL_AMICI_VERSION_STRING";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string getAmiciCommit() const override {
        return "TPL_AMICI_COMMIT_STRING";
    }

    bool hasQuadraticLLH() const override {
        return TPL_QUADRATIC_LLH;
    }

    ObservableScaling getObservableScaling(int iy) const override {
        return observableScalings.at(iy);
    }
};


} // namespace model_TPL_MODELNAME

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
