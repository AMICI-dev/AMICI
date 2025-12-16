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

namespace model_model_nested_events_py {

extern const std::array<std::string_view const, 5> free_parameter_names;
extern const std::array<std::string_view const, 0> fixed_parameter_names;
extern const std::array<std::string_view const, 1> state_names;
extern const std::array<std::string_view const, 1> state_names_solver;
extern const std::array<std::string_view const, 1> observable_names;
extern std::array<const ObservableScaling, 1> observable_scalings;
extern const std::array<std::string_view const, 0> expression_names;
extern const std::array<std::string_view const, 5> free_parameter_ids;
extern const std::array<std::string_view const, 0> fixed_parameter_ids;
extern const std::array<std::string_view const, 1> state_ids;
extern const std::array<std::string_view const, 1> state_ids_solver;
extern const std::array<std::string_view const, 1> observable_ids;
extern const std::array<std::string_view const, 0> expression_ids;
extern std::array<int, 1> state_idxs_solver;

extern void Jy_model_nested_events_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_model_nested_events_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_model_nested_events_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_model_nested_events_py(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_model_nested_events_py(SUNMatrixWrapper &rowvals, int index);






extern void root_model_nested_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);












extern void dxdotdp_explicit_model_nested_events_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_explicit_colptrs_model_nested_events_py(SUNMatrixWrapper &colptrs);
extern void dxdotdp_explicit_rowvals_model_nested_events_py(SUNMatrixWrapper &rowvals);
extern void dxdotdx_explicit_model_nested_events_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdx_explicit_colptrs_model_nested_events_py(SUNMatrixWrapper &colptrs);
extern void dxdotdx_explicit_rowvals_model_nested_events_py(SUNMatrixWrapper &rowvals);
extern void dydx_model_nested_events_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);





extern void sigmay_model_nested_events_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);





extern void x0_model_nested_events_py(realtype *x0, const realtype t, const realtype *p, const realtype *k);

extern void sx0_model_nested_events_py(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);

extern void xdot_model_nested_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_model_nested_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);


extern void stau_model_nested_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie);
extern void deltax_model_nested_events_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old);
extern void deltasx_model_nested_events_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old);
extern void deltaxB_model_nested_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl);
extern void deltaqB_model_nested_events_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB);

extern void x_solver_model_nested_events_py(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_model_nested_events_py(const realtype *p, const realtype *k);


extern std::vector<std::vector<realtype>> explicit_roots_model_nested_events_py(const realtype *p, const realtype *k, const realtype *w);
/**
 * @brief AMICI-generated model subclass.
 */
class Model_model_nested_events_py : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_model_nested_events_py()
        : amici::Model_ODE(
              amici::ModelDimensions{
                  .nx_rdata = 1,
                  .nxtrue_rdata = 1,
                  .nx_solver = 1,
                  .nxtrue_solver = 1,
                  .nx_solver_reinit = 0,
                  .np = 5,
                  .nk = 0,
                  .ny = 1,
                  .nytrue = 1,
                  .nz = 0,
                  .nztrue = 0,
                  .ne = 3,
                  .ne_solver = 2,
                  .nspl = 0,
                  .nw = 0,
                  .ndwdx = 0,
                  .ndwdp = 0,
                  .ndwdw = 0,
                  .ndxdotdw = 0,
                  .ndJydy = std::vector<int>{1},
                  .ndxrdatadxsolver = 0,
                  .ndxrdatadtcl = 0,
                  .ndtotal_cldx_rdata = 0,
                  .nnz = 0,
                  .nJ = 1,
                  .ubw = 1,
                  .lbw = 1,
                  .ndxdotdp_explicit = 2,
                  .ndxdotdx_explicit = 1,
                  .w_recursion_depth = 0,
              },
              amici::SimulationParameters(
                  std::vector<realtype>{}, // fixedParameters
                  std::vector<realtype>{0.10000000000000001, 1000.0, 2.0, 0.80000000000000004, 1.6000000000000001}        // free parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0},   // idlist
              std::vector<int>{},               // z2events
              std::vector<Event>{
                  Event("Heaviside_1", true, true, NAN),
                  Event("Heaviside_2", true, true, NAN),
                  Event("injection", true, true, NAN)
              } // events
          ) {
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_model_nested_events_py(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_model_nested_events_py(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_model_nested_events_py(dJydsigma, iy, p, k, y, sigmay, my);
    }


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {}


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old) override {
        deltax_model_nested_events_py(deltax, t, x, p, k, h, ie, xdot, xdot_old, x_old);
    }


    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old) override {
        deltasx_model_nested_events_py(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau, tcl, x_old);
    }


    void fdeltaxB(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl) override {
        deltaxB_model_nested_events_py(deltaxB, t, x, p, k, h, w, dx, ie, xdot, xdot_old, x_old, xB, tcl);
    }


    void fdeltaqB(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB) override {
        deltaqB_model_nested_events_py(deltaqB, t, x, p, k, h, w, dx, ip, ie, xdot, xdot_old, x_old, xB);
    }


    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_model_nested_events_py(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_model_nested_events_py(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_model_nested_events_py(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_model_nested_events_py(p, k);
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
        dxdotdp_explicit_model_nested_events_py(dxdotdp_explicit, t, x, p, k, h, w);
    }

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdp_explicit_colptrs_model_nested_events_py(colptrs);
    }

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdp_explicit_rowvals_model_nested_events_py(rowvals);
    }


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdx_explicit_model_nested_events_py(dxdotdx_explicit, t, x, p, k, h, w);
    }

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdx_explicit_colptrs_model_nested_events_py(colptrs);
    }

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdx_explicit_rowvals_model_nested_events_py(rowvals);
    }


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dydx_model_nested_events_py(dydx, t, x, p, k, h, w);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {}


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        root_model_nested_events_py(root, t, x, p, k, h, w, tcl);
    }


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_model_nested_events_py(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {
        stau_model_nested_events_py(stau, t, x, p, k, h, w, dx, tcl, sx, ip, ie);
    }

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {
        sx0_model_nested_events_py(sx0, t, x, p, k, ip);
    }

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {}


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {}


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_nested_events_py(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {}


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_model_nested_events_py(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_model_nested_events_py(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_model_nested_events_py(x_solver, x_rdata);
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


    std::vector<std::vector<realtype>> fexplicit_roots(const realtype *p, const realtype *k, const realtype *w) override {
        return explicit_roots_model_nested_events_py(p, k, w);
    }


    std::string get_name() const override {
        return "model_nested_events_py";
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
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string get_amici_version() const override {
        return "1.0.0.dev";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string get_amici_commit() const override {
        return "2b036e909775bc4645c617839fb23d3ed5f07421";
    }

    bool has_quadratic_llh() const override {
        return true;
    }

    ObservableScaling get_observable_scaling(int iy) const override {
        return observable_scalings.at(iy);
    }
};


} // namespace model_model_nested_events_py

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
