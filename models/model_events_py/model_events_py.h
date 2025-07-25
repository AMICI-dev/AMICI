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

namespace model_model_events_py {

extern std::array<const char*, 4> parameterNames;
extern std::array<const char*, 4> fixedParameterNames;
extern std::array<const char*, 3> stateNames;
extern std::array<const char*, 1> observableNames;
extern std::array<const ObservableScaling, 1> observableScalings;
extern std::array<const char*, 1> expressionNames;
extern std::array<const char*, 4> parameterIds;
extern std::array<const char*, 4> fixedParameterIds;
extern std::array<const char*, 3> stateIds;
extern std::array<const char*, 1> observableIds;
extern std::array<const char*, 1> expressionIds;
extern std::array<int, 3> stateIdxsSolver;

extern void Jy_model_events_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydsigma_model_events_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_model_events_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my);
extern void dJydy_colptrs_model_events_py(SUNMatrixWrapper &colptrs, int index);
extern void dJydy_rowvals_model_events_py(SUNMatrixWrapper &rowvals, int index);
extern void Jz_model_events_py(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz);
extern void dJzdsigma_model_events_py(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz);
extern void dJzdz_model_events_py(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz);
extern void Jrz_model_events_py(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz);
extern void dJrzdsigma_model_events_py(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz);
extern void dJrzdz_model_events_py(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz);
extern void root_model_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl);












extern void dxdotdp_explicit_model_events_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_explicit_colptrs_model_events_py(SUNMatrixWrapper &colptrs);
extern void dxdotdp_explicit_rowvals_model_events_py(SUNMatrixWrapper &rowvals);
extern void dxdotdx_explicit_model_events_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdx_explicit_colptrs_model_events_py(SUNMatrixWrapper &colptrs);
extern void dxdotdx_explicit_rowvals_model_events_py(SUNMatrixWrapper &rowvals);
extern void dydx_model_events_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void dydp_model_events_py(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl);
extern void dzdx_model_events_py(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

extern void drzdx_model_events_py(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

extern void sigmay_model_events_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y);
extern void sigmaz_model_events_py(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k);




extern void x0_model_events_py(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void x0_fixedParameters_model_events_py(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs);

extern void sx0_fixedParameters_model_events_py(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs);
extern void xdot_model_events_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_model_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void z_model_events_py(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void rz_model_events_py(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void stau_model_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie);

extern void deltasx_model_events_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old);
extern void deltaxB_model_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl);
extern void deltaqB_model_events_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB);

extern void x_solver_model_events_py(realtype *x_solver, const realtype *x_rdata);












extern std::vector<HermiteSpline> create_splines_model_events_py(const realtype *p, const realtype *k);


extern std::vector<std::vector<realtype>> explicit_roots_model_events_py(const realtype *p, const realtype *k);
/**
 * @brief AMICI-generated model subclass.
 */
class Model_model_events_py : public amici::Model_ODE {
  public:
    /**
     * @brief Default constructor.
     */
    Model_model_events_py()
        : amici::Model_ODE(
              amici::ModelDimensions(
                  3,                            // nx_rdata
                  3,                        // nxtrue_rdata
                  3,                           // nx_solver
                  3,                       // nxtrue_solver
                  0,                    // nx_solver_reinit
                  4,                                  // np
                  4,                                  // nk
                  1,                                  // ny
                  1,                              // nytrue
                  2,                                  // nz
                  2,                              // nztrue
                  6,                              // nevent
                  2,                       // nevent_solver
                  0,                                // nspl
                  1,                          // nobjective
                  1,                                  // nw
                  0,                               // ndwdx
                  0,                               // ndwdp
                  0,                               // ndwdw
                  0,                            // ndxdotdw
                  std::vector<int>{1},                              // ndjydy
                  0,                    // ndxrdatadxsolver
                  0,                        // ndxrdatadtcl
                  0,                        // ndtotal_cldx_rdata
                  0,                                       // nnz
                  3,                                 // ubw
                  3,                                 // lbw
                  true,                                    // pythonGenerated
                  3,                   // ndxdotdp_explicit
                  4,                   // ndxdotdx_explicit
                  0                    // w_recursion_depth
              ),
              amici::SimulationParameters(
                  std::vector<realtype>{4.0, 8.0, 10.0, 4.0}, // fixedParameters
                  std::vector<realtype>{0.5, 2.0, 0.5, 0.5}        // dynamic parameters
              ),
              amici::SecondOrderMode::none,                                  // o2mode
              std::vector<realtype>{1.0, 1.0, 1.0},   // idlist
              std::vector<int>{1, 2},               // z2events
              std::vector<Event>{
                  Event("event_1", false, true, NAN),
                  Event("event_2", false, true, NAN),
                  Event("Heaviside_2", true, true, NAN),
                  Event("Heaviside_3", true, true, NAN),
                  Event("Heaviside_4", true, true, NAN),
                  Event("Heaviside_5", true, true, NAN)
              } // events
          ) {
          }

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    amici::Model *clone() const override {
        return new Model_model_events_py(*this);
    }

    void fJrz(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {
        Jrz_model_events_py(Jrz, iz, p, k, rz, sigmaz);
    }


    void fJy(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        Jy_model_events_py(Jy, iy, p, k, y, sigmay, my);
    }


    void fJz(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {
        Jz_model_events_py(Jz, iz, p, k, z, sigmaz, mz);
    }


    void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {
        dJrzdsigma_model_events_py(dJrzdsigma, iz, p, k, rz, sigmaz);
    }


    void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz) override {
        dJrzdz_model_events_py(dJrzdz, iz, p, k, rz, sigmaz);
    }


    void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydsigma_model_events_py(dJydsigma, iy, p, k, y, sigmay, my);
    }


    void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) override {
        dJzdsigma_model_events_py(dJzdsigma, iz, p, k, z, sigmaz, mz);
    }


    void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz) override {
        dJzdz_model_events_py(dJzdz, iz, p, k, z, sigmaz, mz);
    }


    void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old) override {}


    void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old) override {
        deltasx_model_events_py(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau, tcl, x_old);
    }


    void fdeltaxB(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl) override {
        deltaxB_model_events_py(deltaxB, t, x, p, k, h, dx, ie, xdot, xdot_old, x_old, xB, tcl);
    }


    void fdeltaqB(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB) override {
        deltaqB_model_events_py(deltaqB, t, x, p, k, h, dx, ip, ie, xdot, xdot_old, x_old, xB);
    }


    void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        drzdx_model_events_py(drzdx, ie, t, x, p, k, h);
    }


    void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) override {}


    void fdsigmaydy(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {}


    void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {}


    void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) override {
        dJydy_model_events_py(dJydy, iy, p, k, y, sigmay, my);
    }

    void fdJydy_colptrs(SUNMatrixWrapper &colptrs, int index) override {
        dJydy_colptrs_model_events_py(colptrs, index);
    }

    void fdJydy_rowvals(SUNMatrixWrapper &rowvals, int index) override {
        dJydy_rowvals_model_events_py(rowvals, index);
    }


    std::vector<HermiteSpline> fcreate_splines(const realtype *p, const realtype *k) override {
        return create_splines_model_events_py(p, k);
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
        dxdotdp_explicit_model_events_py(dxdotdp_explicit, t, x, p, k, h, w);
    }

    void fdxdotdp_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdp_explicit_colptrs_model_events_py(colptrs);
    }

    void fdxdotdp_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdp_explicit_rowvals_model_events_py(rowvals);
    }


    void fdxdotdx_explicit(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dxdotdx_explicit_model_events_py(dxdotdx_explicit, t, x, p, k, h, w);
    }

    void fdxdotdx_explicit_colptrs(SUNMatrixWrapper &colptrs) override {
        dxdotdx_explicit_colptrs_model_events_py(colptrs);
    }

    void fdxdotdx_explicit_rowvals(SUNMatrixWrapper &rowvals) override {
        dxdotdx_explicit_rowvals_model_events_py(rowvals);
    }


    void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        dydx_model_events_py(dydx, t, x, p, k, h, w, dwdx);
    }


    void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl) override {
        dydp_model_events_py(dydp, t, x, p, k, h, ip, w, tcl, dtcldp, spl, sspl);
    }


    void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {}


    void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dzdx_model_events_py(dzdx, ie, t, x, p, k, h);
    }


    void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) override {
        root_model_events_py(root, t, x, p, k, h, tcl);
    }


    void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        rz_model_events_py(rz, ie, t, x, p, k, h);
    }


    void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) override {
        sigmay_model_events_py(sigmay, t, p, k, y);
    }


    void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) override {
        sigmaz_model_events_py(sigmaz, t, p, k);
    }


    void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie) override {
        stau_model_events_py(stau, t, x, p, k, h, dx, tcl, sx, ip, ie);
    }

    void fsx0(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) override {}

    void fsx0_fixedParameters(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs) override {
        sx0_fixedParameters_model_events_py(sx0_fixedParameters, t, x0, p, k, ip,  reinitialization_state_idxs);
    }


    void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static) override {}


    void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_events_py(x0, t, p, k);
    }


    void fx0_fixedParameters(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs) override {
        x0_fixedParameters_model_events_py(x0_fixedParameters, t, p, k,  reinitialization_state_idxs);
    }


    void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_model_events_py(xdot, t, x, p, k, h, w);
    }


    void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        y_model_events_py(y, t, x, p, k, h, w);
    }


    void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        z_model_events_py(z, ie, t, x, p, k, h);
    }


    

    void fx_solver(realtype *x_solver, const realtype *x_rdata) override {
        x_solver_model_events_py(x_solver, x_rdata);
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


    std::vector<std::vector<realtype>> fexplicit_roots(const realtype *p, const realtype *k) override {
        return explicit_roots_model_events_py(p, k);
    }


    std::string getName() const override {
        return "model_events_py";
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
        for(auto const idx: stateIdxsSolver) {
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
        return true;
    }

    /**
     * @brief returns the AMICI version that was used to generate the model
     * @return AMICI version string
     */
    std::string getAmiciVersion() const override {
        return "0.33.0";
    }

    /**
     * @brief returns the amici version that was used to generate the model
     * @return AMICI git commit hash
     */
    std::string getAmiciCommit() const override {
        return "d587a622b8295dff051b8cb45d009d9b037f7012";
    }

    bool hasQuadraticLLH() const override {
        return true;
    }

    ObservableScaling getObservableScaling(int iy) const override {
        return observableScalings.at(iy);
    }
};


} // namespace model_model_events_py

} // namespace amici

#endif /* _amici_TPL_MODELNAME_h */
