#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_ode.h"
#include "amici/splinefunctions.h"
#include "amici/event.h"

namespace amici {

class Solver;

namespace model_model_steadystate_py {

extern std::array<const char*, 5> parameter_names;
extern std::array<const char*, 4> fixed_parameter_names;
extern std::array<const char*, 3> state_names;
extern std::array<const char*, 3> observable_names;
extern std::array<const ObservableScaling, 3> observable_scalings;
extern std::array<const char*, 1> expression_names;
extern std::array<const char*, 5> parameter_ids;
extern std::array<const char*, 4> fixed_parameter_ids;
extern std::array<const char*, 3> state_ids;
extern std::array<const char*, 3> observable_ids;
extern std::array<const char*, 1> expression_ids;
extern std::array<int, 3> state_idxs_solver;

extern void Jy_model_steadystate_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_model_steadystate_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_model_steadystate_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_model_steadystate_py(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_model_steadystate_py(SUNMatrixWrapper &rowvals, int index);



















extern void dxdotdp_explicit_model_steadystate_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_explicit_colptrs_model_steadystate_py(SUNMatrixWrapper &colptrs);
extern void dxdotdp_explicit_rowvals_model_steadystate_py(SUNMatrixWrapper &rowvals);
extern void dxdotdx_explicit_model_steadystate_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdx_explicit_colptrs_model_steadystate_py(SUNMatrixWrapper &colptrs);
extern void dxdotdx_explicit_rowvals_model_steadystate_py(SUNMatrixWrapper &rowvals);
extern void dydx_model_steadystate_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);





extern void sigmay_model_steadystate_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);





extern void x0_model_steadystate_py(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_model_steadystate_py(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);

extern void sx0_fixedParameters_model_steadystate_py(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_model_steadystate_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_model_steadystate_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);








extern void x_solver_model_steadystate_py(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_model_steadystate_py(const realtype *p, const realtype *k);



/**
 * @brief AMICI-generated model subclass.
 */
class Model_model_steadystate_py : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_model_steadystate_py()
        : amici::Model_ODE(
              amici::ModelDimensions{
                  .nx_rdata = 3,
                  .nxtrue_rdata = 3,
                  .nx_solver = 3,
                  .nxtrue_solver = 3,
                  .nx_solver_reinit = 0,
                  .np = 5,
                  .nk = 4,
                  .ny = 3,
                  .nytrue = 3,
                  .nz = 0,
                  .nztrue = 0,
                  .ne = 0,
                  .ne_solver = 0,
                  .nspl = 0,
                  .nw = 1,
                  .ndwdx = 0,
                  .ndwdp = 0,
                  .ndwdw = 0,
                  .ndxdotdw = 0,
                  .ndJydy = std::vector<int>{1, 1, 1},
                  .ndxrdatadxsolver = 0,
                  .ndxrdatadtcl = 0,
                  .ndtotal_cldx_rdata = 0,
                  .nnz = 0,
                  .nJ = 1,
                  .ubw = 3,
                  .lbw = 3,
                  .ndxdotdp_explicit = 11,
                  .ndxdotdx_explicit = 9,
                  .w_recursion_depth = 0,
              },
              amici::SimulationParameters(
                  std::vector<realtype>{0.10000000000000001, 0.40000000000000002, 0.69999999999999996, 1.0}, // fixedParameters
                  std::vector<realtype>{1.0, 0.5, 0.40000000000000002, 2.0, 0.10000000000000001}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0},   // idlist
              std::vector<int>{},               // z2events
              std::vector<Event>{} // events
          ) {
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_model_steadystate_py(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_model_steadystate_py(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_model_steadystate_py(dJydsigma, iy, p, k, y, sigmay, my);
    }


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {}


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old) override {}


    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old) override {}


    void fdeltaxB(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl) override {}


    void fdeltaqB(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB) override {}


    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_model_steadystate_py(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_model_steadystate_py(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_model_steadystate_py(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_model_steadystate_py(p, k);
    }

    void fdspline_valuesdp(realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip) override {}

    void fdspline_slopesdp(realtype *dspline_slopesdp, const realtype *p, const realtype *k, const int ip) override {}


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static) override {}

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static) override {}

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static) override {}

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {}

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdp_explicit_model_steadystate_py(dxdotdp_explicit, t, x, p, k, h, w);
    }

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdp_explicit_colptrs_model_steadystate_py(colptrs);
    }

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdp_explicit_rowvals_model_steadystate_py(rowvals);
    }


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdx_explicit_model_steadystate_py(dxdotdx_explicit, t, x, p, k, h, w);
    }

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdx_explicit_colptrs_model_steadystate_py(colptrs);
    }

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdx_explicit_rowvals_model_steadystate_py(rowvals);
    }


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dydx_model_steadystate_py(dydx, t, x, p, k, h, w);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {}


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {}


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_model_steadystate_py(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {}

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {}

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_model_steadystate_py(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);
    }


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {}


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_steadystate_py(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_model_steadystate_py(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_model_steadystate_py(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_model_steadystate_py(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_model_steadystate_py(x_solver, x_rdata);
    }


    void ftotal_cl(realtype *total_cl, const realtype *x_rdata, const realtype *p, const realtype *k) override {}


    void fdx_rdatadx_solver(realtype *dx_rdatadx_solver, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadx_solver_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadx_solver_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdx_rdatadp(realtype *dx_rdatadp, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k, const int ip) override {}


    void fdx_rdatadtcl(realtype *dx_rdatadtcl, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k) override {}

    void fdx_rdatadtcl_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdx_rdatadtcl_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdtotal_cldp(realtype *dtotal_cldp, const realtype *x_rdata, const realtype *p, const realtype *k, const int ip) override {}


    void fdtotal_cldx_rdata(realtype *dtotal_cldx_rdata, const realtype *x_rdata, const realtype *p, const realtype *k, const realtype *tcl) override {}

    void fdtotal_cldx_rdata_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdtotal_cldx_rdata_rowvals(SUNMatrixWrapper &rowvals) override {}


    std::vector<std::vector<realtype>> fexplicit_roots(const realtype *p, const realtype *k) override { return {}; }


    std::string get_name() const override {
        return "model_steadystate_py";
    }

    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    std::vector<std::string> get_parameter_names() const override {
        return std::vector<std::string>(parameter_names.begin(),
                                        parameter_names.end());
    }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    std::vector<std::string> get_state_names() const override {
        return std::vector<std::string>(state_names.begin(), state_names.end());
    }

    /**
     * @brief Get names of the solver states
     * @return the names
     */
    std::vector<std::string> get_state_names_solver() const override {
        std::vector<std::string> result;
        result.reserve(state_idxs_solver.size());
        for(auto const idx: state_idxs_solver) {
            result.push_back(state_names[idx]);
        }
        return result;
    }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    std::vector<std::string> get_fixed_parameter_names() const override {
        return std::vector<std::string>(fixed_parameter_names.begin(),
                                        fixed_parameter_names.end());
    }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    std::vector<std::string> get_observable_names() const override {
        return std::vector<std::string>(observable_names.begin(),
                                        observable_names.end());
    }

    /**
     * @brief Get names of model expressions
     * @return Expression names
     */
    std::vector<std::string> get_expression_names() const override {
        return std::vector<std::string>(expression_names.begin(),
                                        expression_names.end());
    }

    /**
     * @brief Get ids of the model parameters
     * @return the ids
     */
    std::vector<std::string> get_parameter_ids() const override {
        return std::vector<std::string>(parameter_ids.begin(),
                                        parameter_ids.end());
    }

    /**
     * @brief Get ids of the model states
     * @return the ids
     */
    std::vector<std::string> get_state_ids() const override {
        return std::vector<std::string>(state_ids.begin(), state_ids.end());
    }

    /**
     * @brief Get ids of the solver states
     * @return the ids
     */
    std::vector<std::string> get_state_ids_solver() const override {
        std::vector<std::string> result;
        result.reserve(state_idxs_solver.size());
        for(auto idx: state_idxs_solver) {
            result.push_back(state_ids[idx]);
        }
        return result;
    }

    /**
     * @brief Get ids of the fixed model parameters
     * @return the ids
     */
    std::vector<std::string> get_fixed_parameter_ids() const override {
        return std::vector<std::string>(fixed_parameter_ids.begin(),
                                        fixed_parameter_ids.end());
    }

    /**
     * @brief Get ids of the observables
     * @return the ids
     */
    std::vector<std::string> get_observable_ids() const override {
        return std::vector<std::string>(observable_ids.begin(),
                                        observable_ids.end());
    }

    /**
     * @brief Get IDs of model expressions
     * @return Expression IDs
     */
    std::vector<std::string> get_expression_ids() const override {
        return std::vector<std::string>(expression_ids.begin(),
                                        expression_ids.end());
    }

    /**
     * @brief function indicating whether reinitialization of states depending
     * on fixed parameters is permissible
     * @return flag indicating whether reinitialization of states depending on
     * fixed parameters is permissible
     */
    bool is_fixed_parameter_state_reinitialization_allowed() const override {
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string get_amici_version() const override {
        return "0.34.1";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string get_amici_commit() const override {
        return "f005fac9e2de7c3c90be2ac55d4ad165471ed1e7";
    }

    bool has_quadratic_llh() const override {
        return true;
    }

    ObservableScaling get_observable_scaling(int iy) const override {
        return observable_scalings.at(iy);
    }
};


} // namespace model_model_steadystate_py

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
