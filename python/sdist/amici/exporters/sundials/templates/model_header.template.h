#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_TPL_MODEL_TYPE_LOWER.h"
#include "amici/splinefunctions.h"
#include "amici/event.h"

namespace amici {

class Solver;

namespace model_TPL_MODELNAME {

extern const std::array<std::string_view const, TPL_NP> free_parameter_names;
extern const std::array<std::string_view const, TPL_NK> fixed_parameter_names;
extern const std::array<std::string_view const, TPL_NX_RDATA> state_names;
extern const std::array<std::string_view const, TPL_NX_SOLVER> state_names_solver;
extern const std::array<std::string_view const, TPL_NY> observable_names;
extern std::array<const ObservableScaling, TPL_NY> observable_scalings;
extern const std::array<std::string_view const, TPL_NW> expression_names;
extern const std::array<std::string_view const, TPL_NP> free_parameter_ids;
extern const std::array<std::string_view const, TPL_NK> fixed_parameter_ids;
extern const std::array<std::string_view const, TPL_NX_RDATA> state_ids;
extern const std::array<std::string_view const, TPL_NX_SOLVER> state_ids_solver;
extern const std::array<std::string_view const, TPL_NY> observable_ids;
extern const std::array<std::string_view const, TPL_NW> expression_ids;
extern std::array<int, TPL_NX_SOLVER> state_idxs_solver;

TPL_JY_DEF
TPL_DJYDSIGMA_DEF
TPL_DJYDY_DEF
TPL_DJYDY_COLPTRS_DEF
TPL_DJYDY_ROWVALS_DEF
TPL_JZ_DEF
TPL_DJZDSIGMA_DEF
TPL_DJZDZ_DEF
TPL_JRZ_DEF
TPL_DJRZDSIGMA_DEF
TPL_DJRZDZ_DEF
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
TPL_DZDX_DEF
TPL_DZDP_DEF
TPL_DRZDX_DEF
TPL_DRZDP_DEF
TPL_SIGMAY_DEF
TPL_SIGMAZ_DEF
TPL_DSIGMAYDP_DEF
TPL_DSIGMAYDY_DEF
TPL_DSIGMAZDP_DEF
TPL_W_DEF
TPL_X0_DEF
TPL_X0_FIXEDPARAMETERS_DEF
TPL_SX0_DEF
TPL_SX0_FIXEDPARAMETERS_DEF
TPL_XDOT_DEF
TPL_Y_DEF
TPL_Z_DEF
TPL_RZ_DEF
TPL_STAU_DEF
TPL_DELTAX_DEF
TPL_DELTASX_DEF
TPL_DELTAXB_DEF
TPL_DELTAQB_DEF
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
TPL_CREATE_SPLINES_DEF
TPL_DSPLINE_VALUESDP_DEF
TPL_DSPLINE_SLOPESDP_DEF
TPL_EXPLICIT_ROOTS_DEF
/**
 * @brief AMICI-generated model subclass.
 */
class Model_TPL_MODELNAME : public amici::Model_TPL_MODEL_TYPE_UPPER {
  public:
    /**
     * @brief Default constructor.
     */
    Model_TPL_MODELNAME()
        : amici::Model_TPL_MODEL_TYPE_UPPER(
              amici::ModelDimensions{
                  .nx_rdata = TPL_NX_RDATA,
                  .nxtrue_rdata = TPL_NXTRUE_RDATA,
                  .nx_solver = TPL_NX_SOLVER,
                  .nxtrue_solver = TPL_NXTRUE_SOLVER,
                  .nx_solver_reinit = TPL_NX_SOLVER_REINIT,
                  .np = TPL_NP,
                  .nk = TPL_NK,
                  .ny = TPL_NY,
                  .nytrue = TPL_NYTRUE,
                  .nz = TPL_NZ,
                  .nztrue = TPL_NZTRUE,
                  .ne = TPL_NEVENT,
                  .ne_solver = TPL_NEVENT_SOLVER,
                  .nspl = TPL_NSPL,
                  .nw = TPL_NW,
                  .ndwdx = TPL_NDWDX,
                  .ndwdp = TPL_NDWDP,
                  .ndwdw = TPL_NDWDW,
                  .ndxdotdw = TPL_NDXDOTDW,
                  .ndJydy = TPL_NDJYDY,
                  .ndxrdatadxsolver = TPL_NDXRDATADXSOLVER,
                  .ndxrdatadtcl = TPL_NDXRDATADTCL,
                  .ndtotal_cldx_rdata = TPL_NDTOTALCLDXRDATA,
                  .nnz = 0,
                  .nJ = TPL_NOBJECTIVE,
                  .ubw = TPL_UBW,
                  .lbw = TPL_LBW,
                  .ndxdotdp_explicit = TPL_NDXDOTDP_EXPLICIT,
                  .ndxdotdx_explicit = TPL_NDXDOTDX_EXPLICIT,
                  .w_recursion_depth = TPL_W_RECURSION_DEPTH,
              },
              amici::SimulationParameters(
                  std::vector<realtype>{TPL_FIXED_PARAMETERS}, // fixedParameters
                  std::vector<realtype>{TPL_FREE_PARAMETERS}        // free parameters
              ),
              TPL_O2MODE,                                  // o2mode
              std::vector<realtype>{TPL_ID},   // idlist
              std::vector<int>{TPL_Z2EVENT},               // z2events
              std::vector<Event>{TPL_EVENT_LIST_INITIALIZER} // events
          ) {
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_TPL_MODELNAME(*this);
    }

    TPL_JRZ_IMPL

    TPL_JY_IMPL

    TPL_JZ_IMPL

    TPL_DJRZDSIGMA_IMPL

    TPL_DJRZDZ_IMPL

    TPL_DJYDSIGMA_IMPL

    TPL_DJZDSIGMA_IMPL

    TPL_DJZDZ_IMPL

    TPL_DELTAX_IMPL

    TPL_DELTASX_IMPL

    TPL_DELTAXB_IMPL

    TPL_DELTAQB_IMPL

    TPL_DRZDP_IMPL

    TPL_DRZDX_IMPL

    TPL_DSIGMAYDP_IMPL

    TPL_DSIGMAYDY_IMPL

    TPL_DSIGMAZDP_IMPL

    TPL_DJYDY_IMPL
    TPL_DJYDY_COLPTRS_IMPL
    TPL_DJYDY_ROWVALS_IMPL

    TPL_CREATE_SPLINES_IMPL
    TPL_DSPLINE_VALUESDP_IMPL
    TPL_DSPLINE_SLOPESDP_IMPL

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

    TPL_DZDP_IMPL

    TPL_DZDX_IMPL

    TPL_ROOT_IMPL

    TPL_RZ_IMPL

    TPL_SIGMAY_IMPL

    TPL_SIGMAZ_IMPL

    TPL_STAU_IMPL
    TPL_SX0_IMPL
    TPL_SX0_FIXEDPARAMETERS_IMPL

    TPL_W_IMPL

    TPL_X0_IMPL

    TPL_X0_FIXEDPARAMETERS_IMPL

    TPL_XDOT_IMPL

    TPL_Y_IMPL

    TPL_Z_IMPL

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

    TPL_EXPLICIT_ROOTS_IMPL

    std::string get_name() const override {
        return "TPL_MODELNAME";
    }

    /**
     * @brief Get names of the free model parameters
     * @return the names
     */
    std::span<std::string_view const> get_free_parameter_names() const override {
        return free_parameter_names;
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    std::span<std::string_view const> get_state_names() const override {
        return state_names;
    }

    /**
     * @brief Get names of the solver states
     * @return the names
     */
    std::span<std::string_view const> get_state_names_solver() const override {
        return state_names_solver;
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    std::span<std::string_view const> get_fixed_parameter_names() const override {
        return fixed_parameter_names;
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    std::span<std::string_view const> get_observable_names() const override {
        return observable_names;
    }

    /**
     * @brief Get names of model expressions
     * @return Expression names
     */
    std::span<std::string_view const> get_expression_names() const override {
        return expression_names;
    }

    /**
     * @brief Get ids of the free model parameters
     * @return the ids
     */
    std::span<std::string_view const> get_free_parameter_ids() const override {
        return free_parameter_ids;
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    std::span<std::string_view const> get_state_ids() const override {
        return state_ids;
    }

    /**
     * @brief Get ids of the solver states
     * @return the ids
     */
    std::span<std::string_view const> get_state_ids_solver() const override {
        return state_ids_solver;
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    std::span<std::string_view const> get_fixed_parameter_ids() const override {
        return fixed_parameter_ids;
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    std::span<std::string_view const> get_observable_ids() const override {
        return observable_ids;
    }

    /**
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    std::span<std::string_view const> get_expression_ids() const override {
        return expression_ids;
    }

    /**
     * @brief function indicating whether reinitialization of states depending
     * on fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    bool is_fixed_parameter_state_reinitialization_allowed() const override {
        return TPL_REINIT_FIXPAR_INITCOND;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string get_amici_version() const override {
        return "TPL_AMICI_VERSION_STRING";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string get_amici_commit() const override {
        return "TPL_AMICI_COMMIT_STRING";
    }

    bool has_quadratic_llh() const override {
        return TPL_QUADRATIC_LLH;
    }

    ObservableScaling get_observable_scaling(int iy) const override {
        return observable_scalings.at(iy);
    }
};


} // namespace model_TPL_MODELNAME

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
