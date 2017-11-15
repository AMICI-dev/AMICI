#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include <math.h>
#include <include/amici_solver_cvodes.h>
#include <include/amici_model_ode.h>

namespace amici {
class UserData;
class Solver;
}


#define pi M_PI

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

amici::UserData getUserData();
amici::Solver *getSolver();
amici::Model *getModel();

class Model_model_events : public amici::Model_ODE {
public:
    Model_model_events() : amici::Model_ODE(4,
                    3,
                    3,
                    4,
                    1,
                    1,
                    2,
                    2,
                    4,
                    1,
                    0,
                    0,
                    0,
                    4,
                    0,
                    1,
                    amici::AMICI_O2MODE_NONE)
    {
        z2event = new int[2] {1, 2,};
        idlist = new realtype[3] {0, 0, 0,};
    }

    void J_model_events(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx);
    void model_J = &J_model_events;

    void JB_model_events(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
    void model_JB = &JB_model_events;

    void JDiag_model_events(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
    void model_JDiag = &JDiag_model_events;

    void JSparse_model_events(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
    void model_JSparse = &JSparse_model_events;

    void JSparseB_model_events(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
    void model_JSparseB = &JSparseB_model_events;

    void Jrz_model_events(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
    void model_Jrz = &Jrz_model_events;

    void Jv_model_events(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx);
    void model_Jv = &Jv_model_events;

    void JvB_model_events(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *vB, const realtype *w, const realtype *dwdx);
    void model_JvB = &JvB_model_events;

    void Jy_model_events(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    void model_Jy = &Jy_model_events;

    void Jz_model_events(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    void model_Jz = &Jz_model_events;

    void dJrzdsigma_model_events(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
    void model_dJrzdsigma = &dJrzdsigma_model_events;

    void dJrzdz_model_events(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
    void model_dJrzdz = &dJrzdz_model_events;

    void dJydsigma_model_events(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    void model_dJydsigma = &dJydsigma_model_events;

    void dJydy_model_events(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    void model_dJydy = &dJydy_model_events;

    void dJzdsigma_model_events(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    void model_dJzdsigma = &dJzdsigma_model_events;

    void dJzdz_model_events(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    void model_dJzdz = &dJzdz_model_events;

    void deltaqB_model_events(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB, const realtype *qBdot);
    void model_deltaqB = &deltaqB_model_events;

    void deltasx_model_events(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau);
    void model_deltasx = &deltasx_model_events;

    void deltax_model_events(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old);
    void model_deltax = &deltax_model_events;

    void deltaxB_model_events(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
    void model_deltaxB = &deltaxB_model_events;

    void drzdp_model_events(double *drzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
    void model_drzdp = &drzdp_model_events;

    void drzdx_model_events(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_drzdx = &drzdx_model_events;

    void dsigma_ydp_model_events(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip);
    void model_dsigma_ydp = &dsigma_ydp_model_events;

    void dsigma_zdp_model_events(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip);
    void model_dsigma_zdp = &dsigma_zdp_model_events;

    void dwdp_model_events(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
    void model_dwdp = &dwdp_model_events;

    void dwdx_model_events(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
    void model_dwdx = &dwdx_model_events;

    void dxdotdp_model_events(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp);
    void model_dxdotdp = &dxdotdp_model_events;

    void dydp_model_events(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
    void model_dydp = &dydp_model_events;

    void dydx_model_events(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_dydx = &dydx_model_events;

    void dzdp_model_events(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
    void model_dzdp = &dzdp_model_events;

    void dzdx_model_events(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_dzdx = &dzdx_model_events;

    void qBdot_model_events(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp);
    void model_qBdot = &qBdot_model_events;

    void root_model_events(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_root = &root_model_events;

    void rz_model_events(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_rz = &rz_model_events;

    void sigma_y_model_events(double *sigmay, const realtype t, const realtype *p, const realtype *k);
    void model_sigma_y = &sigma_y_model_events;

    void sigma_z_model_events(double *sigmaz, const realtype t, const realtype *p, const realtype *k);
    void model_sigma_z = &sigma_z_model_events;

    void srz_model_events(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
    void model_srz = &srz_model_events;

    void stau_model_events(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie);
    void model_stau = &stau_model_events;

    void sx0_model_events(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip);
    void model_sx0 = &sx0_model_events;

    void sxdot_model_events(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp);
    void model_sxdot = &sxdot_model_events;

    void sz_model_events(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
    void model_sz = &sz_model_events;

    void w_model_events(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_w = &w_model_events;

    void x0_model_events(realtype *x0, const realtype t, const realtype *p, const realtype *k);
    void model_x0 = &x0_model_events;

    void xBdot_model_events(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
    void model_xBdot = &xBdot_model_events;

    void xdot_model_events(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
    void model_xdot = &xdot_model_events;

    void y_model_events(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_y = &y_model_events;

    void z_model_events(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    void model_z = &z_model_events;

    amici::Solver *getSolver() override {
        return new amici::CVodeSolver();
    }

#endif /* _amici_wrapfunctions_h */
