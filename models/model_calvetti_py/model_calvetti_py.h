#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include <gsl/gsl-lite.hpp>

#include "amici/model_dae.h"
#include "amici/splinefunctions.h"
#include "amici/event.h"

namespace amici {

class Solver;

namespace model_model_calvetti_py {

extern const std::array<std::string_view const, 0> free_parameter_names;
extern const std::array<std::string_view const, 6> fixed_parameter_names;
extern const std::array<std::string_view const, 6> state_names;
extern const std::array<std::string_view const, 6> state_names_solver;
extern const std::array<std::string_view const, 6> observable_names;
extern std::array<const ObservableScaling, 6> observable_scalings;
extern const std::array<std::string_view const, 16> expression_names;
extern const std::array<std::string_view const, 0> free_parameter_ids;
extern const std::array<std::string_view const, 6> fixed_parameter_ids;
extern const std::array<std::string_view const, 6> state_ids;
extern const std::array<std::string_view const, 6> state_ids_solver;
extern const std::array<std::string_view const, 6> observable_ids;
extern const std::array<std::string_view const, 16> expression_ids;
extern std::array<int, 6> state_idxs_solver;

extern void Jy_model_calvetti_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_model_calvetti_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_model_calvetti_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_model_calvetti_py(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_model_calvetti_py(SUNMatrixWrapper &rowvals, int index);






extern void root_model_calvetti_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl);



extern void dwdx_model_calvetti_py(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static);
extern void dwdx_colptrs_model_calvetti_py(SUNMatrixWrapper &colptrs);
extern void dwdx_rowvals_model_calvetti_py(SUNMatrixWrapper &rowvals);
extern void dwdw_model_calvetti_py(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static);
extern void dwdw_colptrs_model_calvetti_py(SUNMatrixWrapper &colptrs);
extern void dwdw_rowvals_model_calvetti_py(SUNMatrixWrapper &rowvals);
extern void dxdotdw_model_calvetti_py(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w);
extern void dxdotdw_colptrs_model_calvetti_py(SUNMatrixWrapper &colptrs);
extern void dxdotdw_rowvals_model_calvetti_py(SUNMatrixWrapper &rowvals);



extern void dxdotdx_explicit_model_calvetti_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w);
extern void dxdotdx_explicit_colptrs_model_calvetti_py(SUNMatrixWrapper &colptrs);
extern void dxdotdx_explicit_rowvals_model_calvetti_py(SUNMatrixWrapper &rowvals);
extern void dydx_model_calvetti_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);





extern void sigmay_model_calvetti_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);




extern void w_model_calvetti_py(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static);
extern void x0_model_calvetti_py(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_model_calvetti_py(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);


extern void xdot_model_calvetti_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w);
extern void y_model_calvetti_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);








extern void x_solver_model_calvetti_py(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_model_calvetti_py(const realtype *p, const realtype *k);


extern std::vector<std::vector<realtype>> explicit_roots_model_calvetti_py(const realtype *p, const realtype *k, const realtype *w);
/**
 * @brief AMICI-generated model subclass.
 */
class Model_model_calvetti_py : public amici::Model_DAE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_model_calvetti_py()
        : amici::Model_DAE(
              amici::ModelDimensions{
                  .nx_rdata = 6,
                  .nxtrue_rdata = 6,
                  .nx_solver = 6,
                  .nxtrue_solver = 6,
                  .nx_solver_reinit = 0,
                  .np = 0,
                  .nk = 6,
                  .ny = 6,
                  .nytrue = 6,
                  .nz = 0,
                  .nztrue = 0,
                  .ne = 4,
                  .ne_solver = 0,
                  .nspl = 0,
                  .nw = 16,
                  .ndwdx = 15,
                  .ndwdp = 0,
                  .ndwdw = 16,
                  .ndxdotdw = 7,
                  .ndJydy = std::vector<int>{1, 1, 1, 1, 1, 1},
                  .ndxrdatadxsolver = 0,
                  .ndxrdatadtcl = 0,
                  .ndtotal_cldx_rdata = 0,
                  .nnz = 0,
                  .nJ = 1,
                  .ubw = 6,
                  .lbw = 6,
                  .ndxdotdp_explicit = 0,
                  .ndxdotdx_explicit = 5,
                  .w_recursion_depth = 2,
              },
              amici::SimulationParameters(
                  std::vector<realtype>{0.28999999999999998, 0.73999999999999999, 0.44, 0.080000000000000002, 0.27000000000000002, 0.17999999999999999}, // fixedParameters
                  std::vector<realtype>{}        // free parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0, 0.0, 0.0, 0.0},   // idlist
              std::vector<int>{},               // z2events
              std::vector<Event>{
                  Event("Heaviside_0", true, true, NAN),
                  Event("Heaviside_1", true, true, NAN),
                  Event("Heaviside_2", true, true, NAN),
                  Event("Heaviside_3", true, true, NAN)
              } // events
          ) {
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_model_calvetti_py(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_model_calvetti_py(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {}


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_model_calvetti_py(dJydsigma, iy, p, k, y, sigmay, my);
    }


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {}


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {}


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old) override {}


    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old) override {}


    void fdeltaxB(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl) override {}


    void fdeltaqB(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB) override {}


    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_model_calvetti_py(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_model_calvetti_py(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_model_calvetti_py(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_model_calvetti_py(p, k);
    }

    void fdspline_valuesdp(realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip) override {}

    void fdspline_slopesdp(realtype *dspline_slopesdp, const realtype *p, const realtype *k, const int ip) override {}


    void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static) override {}

    void fdwdp_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdwdp_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static) override {
        dwdx_model_calvetti_py(dwdx, t, x, p, k, h, w, tcl, spl, include_static);
    }

    void fdwdx_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdx_colptrs_model_calvetti_py(colptrs);
    }

    void fdwdx_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdx_rowvals_model_calvetti_py(rowvals);
    }


    void fdwdw(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static) override {
        dwdw_model_calvetti_py(dwdw, t, x, p, k, h, w, tcl, include_static);
    }

    void fdwdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dwdw_colptrs_model_calvetti_py(colptrs);
    }

    void fdwdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dwdw_rowvals_model_calvetti_py(rowvals);
    }


    void fdxdotdw(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) override {
        dxdotdw_model_calvetti_py(dxdotdw, t, x, p, k, h, dx, w);
    }

    void fdxdotdw_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdw_colptrs_model_calvetti_py(colptrs);
    }

    void fdxdotdw_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdw_rowvals_model_calvetti_py(rowvals);
    }


    void fdxdotdp_explicit(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) override {}

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {}

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {}


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) override {
        dxdotdx_explicit_model_calvetti_py(dxdotdx_explicit, t, x, p, k, h, dx, w);
    }

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdx_explicit_colptrs_model_calvetti_py(colptrs);
    }

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdx_explicit_rowvals_model_calvetti_py(rowvals);
    }


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dydx_model_calvetti_py(dydx, t, x, p, k, h, w);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {}


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) override {
        root_model_calvetti_py(root, t, x, p, k, h, w, tcl);
    }


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_model_calvetti_py(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {}


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {}

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {}

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {}


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {
        w_model_calvetti_py(w, t, x, p, k, h, tcl, spl, include_static);
    }


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_calvetti_py(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_model_calvetti_py(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) override {
        xdot_model_calvetti_py(xdot, t, x, p, k, h, dx, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_model_calvetti_py(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {}


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_model_calvetti_py(x_solver, x_rdata);
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
        return explicit_roots_model_calvetti_py(p, k, w);
    }


    std::string get_name() const override {
        return "model_calvetti_py";
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


} // namespace model_model_calvetti_py

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
