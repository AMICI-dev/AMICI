#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include <math.h>
#include <include/amici_defines.h>
#include <sundials/sundials_sparse.h> //SlsMat definition
#include <include/udata.h>
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

amici::Model *getModel(const amici::UserData *udata);
void getModelDims(int *nx, int *nk, int *np);

extern void J_model_steadystate(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JB_model_steadystate(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JDiag_model_steadystate(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_model_steadystate(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparseB_model_steadystate(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void Jrz_model_steadystate(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void Jv_model_steadystate(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx);
extern void JvB_model_steadystate(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx);
extern void Jy_model_steadystate(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void Jz_model_steadystate(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void dJrzdsigma_model_steadystate(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void dJrzdz_model_steadystate(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void dJydsigma_model_steadystate(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydy_model_steadystate(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJzdsigma_model_steadystate(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void dJzdz_model_steadystate(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void deltaqB_model_steadystate(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
extern void deltasx_model_steadystate(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau);
extern void deltax_model_steadystate(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old);
extern void deltaxB_model_steadystate(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
extern void drzdp_model_steadystate(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void drzdx_model_steadystate(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dsigma_ydp_model_steadystate(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip);
extern void dsigma_zdp_model_steadystate(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip);
extern void dwdp_model_steadystate(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dwdx_model_steadystate(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_model_steadystate(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp);
extern void dydp_model_steadystate(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void dydx_model_steadystate(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dzdp_model_steadystate(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void dzdx_model_steadystate(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void qBdot_model_steadystate(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp);
extern void root_model_steadystate(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void rz_model_steadystate(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void sigma_y_model_steadystate(double *sigmay, const realtype t, const realtype *p, const realtype *k);
extern void sigma_z_model_steadystate(double *sigmaz, const realtype t, const realtype *p, const realtype *k);
extern void srz_model_steadystate(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
extern void stau_model_steadystate(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie);
extern void sx0_model_steadystate(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip);
extern void sxdot_model_steadystate(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp);
extern void sz_model_steadystate(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
extern void w_model_steadystate(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void x0_model_steadystate(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void xBdot_model_steadystate(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void xdot_model_steadystate(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_model_steadystate(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void z_model_steadystate(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

class Model_model_steadystate : public amici::Model_ODE {
public:
    Model_model_steadystate(const amici::UserData *udata) : amici::Model_ODE(5,
                    3,
                    3,
                    4,
                    3,
                    3,
                    0,
                    0,
                    0,
                    1,
                    2,
                    2,
                    1,
                    9,
                    2,
                    2,
                    amici::AMICI_O2MODE_NONE,
                    std::vector<realtype>(udata->unp(),udata->unp()+5),
                    std::vector<realtype>(udata->k(),udata->k()+4),
                    udata->plist(),
                    std::vector<realtype>{0, 0, 0},
                    std::vector<int>{})
                    {};

    virtual void model_J(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        J_model_steadystate(J, t, x, p, k, h, w, dwdx);
    }

    virtual int model_JB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JB_model_steadystate(JB, t, x, p, k, h, xB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JDiag(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JDiag_model_steadystate(JDiag, t, x, p, k, h, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JSparse(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_model_steadystate(JSparse, t, x, p, k, h, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JSparseB(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_model_steadystate(JSparseB, t, x, p, k, h, xB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_Jrz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        Jrz_model_steadystate(nllh, iz, p, k, rz, sigmaz);
    }

    virtual int model_Jv(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) override {
        Jv_model_steadystate(Jv, t, x, p, k, h, v, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JvB(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) override {
        JvB_model_steadystate(JvB, t, x, p, k, h, xB, vB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_Jy(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        Jy_model_steadystate(nllh, iy, p, k, y, sigmay, my);
    }

    virtual void model_Jz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        Jz_model_steadystate(nllh, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_dJrzdsigma(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        dJrzdsigma_model_steadystate(dJrzdsigma, iz, p, k, rz, sigmaz);
    }

    virtual void model_dJrzdz(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        dJrzdz_model_steadystate(dJrzdz, iz, p, k, rz, sigmaz);
    }

    virtual void model_dJydsigma(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydsigma_model_steadystate(dJydsigma, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJydy(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydy_model_steadystate(dJydy, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJzdsigma(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        dJzdsigma_model_steadystate(dJzdsigma, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_dJzdz(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        dJzdz_model_steadystate(dJzdz, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_deltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
        deltaqB_model_steadystate(deltaqB, t, x, p, k, h, ip, ie, xdot, xdot_old, xB);
    }

    virtual void model_deltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) override {
        deltasx_model_steadystate(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau);
    }

    virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {
        deltax_model_steadystate(deltax, t, x, p, k, h, ie, xdot, xdot_old);
    }

    virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
        deltaxB_model_steadystate(deltaxB, t, x, p, k, h, ie, xdot, xdot_old, xB);
    }

    virtual void model_drzdp(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        drzdp_model_steadystate(drzdp, ie, t, x, p, k, h, ip);
    }

    virtual void model_drzdx(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        drzdx_model_steadystate(drzdx, ie, t, x, p, k, h);
    }

    virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
        dsigma_ydp_model_steadystate(dsigmaydp, t, p, k, ip);
    }

    virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
        dsigma_zdp_model_steadystate(dsigmazdp, t, p, k, ip);
    }

    virtual void model_dwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdp_model_steadystate(dwdp, t, x, p, k, h, w);
    }

    virtual void model_dwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdx_model_steadystate(dwdx, t, x, p, k, h, w);
    }

    virtual int model_dxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) override {
        dxdotdp_model_steadystate(dxdotdp, t, x, p, k, h, ip, w, dwdp);
        return AMICI_SUCCESS;
    }

    virtual void model_dydp(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        dydp_model_steadystate(dydp, t, x, p, k, h, ip);
    }

    virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dydx_model_steadystate(dydx, t, x, p, k, h);
    }

    virtual void model_dzdp(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        dzdp_model_steadystate(dzdp, ie, t, x, p, k, h, ip);
    }

    virtual void model_dzdx(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dzdx_model_steadystate(dzdx, ie, t, x, p, k, h);
    }

    virtual int model_qBdot(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) override {
        qBdot_model_steadystate(qBdot, ip, t, x, p, k, h, xB, w, dwdp);
        return AMICI_SUCCESS;
    }

    virtual int model_root(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        root_model_steadystate(root, t, x, p, k, h);
        return AMICI_SUCCESS;
    }

    virtual void model_rz(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        rz_model_steadystate(rz, ie, t, x, p, k, h);
    }

    virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) override {
        sigma_y_model_steadystate(sigmay, t, p, k);
    }

    virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) override {
        sigma_z_model_steadystate(sigmaz, t, p, k);
    }

    virtual void model_srz(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
        srz_model_steadystate(srz, ie, t, x, p, k, h, sx, ip);
    }

    virtual void model_stau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) override {
        stau_model_steadystate(stau, t, x, p, k, h, sx, ip, ie);
    }

    virtual void model_sx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) override {
        sx0_model_steadystate(sx0, t, x0, p, k, ip);
    }

    virtual int model_sxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) override {
        sxdot_model_steadystate(sxdot, t, x, p, k, h, ip, sx, w, dwdx, J, dxdotdp);
        return AMICI_SUCCESS;
    }

    virtual void model_sz(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
        sz_model_steadystate(sz, ie, t, x, p, k, h, sx, ip);
    }

    virtual void model_w(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        w_model_steadystate(w, t, x, p, k, h);
    }

    virtual void model_x0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_steadystate(x0, t, p, k);
    }

    virtual int model_xBdot(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        xBdot_model_steadystate(xBdot, t, x, p, k, h, xB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_xdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_model_steadystate(xdot, t, x, p, k, h, w);
    }

    virtual void model_y(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        y_model_steadystate(y, t, x, p, k, h);
    }

    virtual void model_z(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        z_model_steadystate(z, ie, t, x, p, k, h);
    }

};

#endif /* _amici_wrapfunctions_h */
