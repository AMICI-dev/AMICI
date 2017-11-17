#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include <math.h>
#include <include/amici_defines.h>
#include <include/udata.h>
#include <include/amici_solver_idas.h>
#include <include/amici_model_dae.h>

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

extern void J_model_robertson(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx);
extern void JB_model_robertson(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx);
extern void JDiag_model_robertson(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx);
extern void JSparse_model_robertson(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx);
extern void JSparseB_model_robertson(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx);
extern void Jrz_model_robertson(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void Jv_model_robertson(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *v, const realtype *w, const realtype *dwdx);
extern void JvB_model_robertson(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *vB, const realtype *w, const realtype *dwdx);
extern void Jy_model_robertson(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void Jz_model_robertson(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void M_model_robertson(realtype *M, const realtype t, const realtype *x, const realtype *p, const realtype *k);
extern void dJrzdsigma_model_robertson(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void dJrzdz_model_robertson(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
extern void dJydsigma_model_robertson(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydy_model_robertson(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJzdsigma_model_robertson(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void dJzdz_model_robertson(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
extern void deltaqB_model_robertson(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
extern void deltasx_model_robertson(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau);
extern void deltax_model_robertson(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old);
extern void deltaxB_model_robertson(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
extern void drzdp_model_robertson(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void drzdx_model_robertson(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dsigma_ydp_model_robertson(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip);
extern void dsigma_zdp_model_robertson(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip);
extern void dwdp_model_robertson(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dwdx_model_robertson(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_model_robertson(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *dx, const realtype *w, const realtype *dwdp);
extern void dydp_model_robertson(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void dydx_model_robertson(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dzdp_model_robertson(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void dzdx_model_robertson(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void qBdot_model_robertson(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdp);
extern void root_model_robertson(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx);
extern void rz_model_robertson(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void sigma_y_model_robertson(double *sigmay, const realtype t, const realtype *p, const realtype *k);
extern void sigma_z_model_robertson(double *sigmaz, const realtype t, const realtype *p, const realtype *k);
extern void srz_model_robertson(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
extern void stau_model_robertson(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie);
extern void sx0_model_robertson(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip);
extern void sxdot_model_robertson(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *dx, const realtype *sx, const realtype *sdx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *M, const realtype *dxdotdp);
extern void sz_model_robertson(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip);
extern void w_model_robertson(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void x0_model_robertson(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void xBdot_model_robertson(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx);
extern void xdot_model_robertson(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w);
extern void y_model_robertson(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void z_model_robertson(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

class Model_model_robertson : public amici::Model_DAE {
public:
    Model_model_robertson(const amici::UserData *udata) : amici::Model_DAE(3,
                    3,
                    3,
                    1,
                    3,
                    3,
                    0,
                    0,
                    0,
                    1,
                    1,
                    2,
                    1,
                    9,
                    2,
                    2,
                    amici::AMICI_O2MODE_NONE,
                    std::vector<realtype>(udata->unp(),udata->unp()+3),
                    std::vector<realtype>(udata->k(),udata->k()+1),
                    udata->plist(),
                    std::vector<realtype>{1, 1, 0},
                    std::vector<int>{})
                    {};

    virtual void model_J(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) override {
        J_model_robertson(J, t, x, p, k, h, cj, dx, w, dwdx);
    }

    virtual int model_JB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) override {
        JB_model_robertson(JB, t, x, p, k, h, cj, xB, dx, dxB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JDiag(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) override {
        JDiag_model_robertson(JDiag, t, x, p, k, h, cj, dx, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JSparse(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) override {
        JSparse_model_robertson(J, t, x, p, k, h, cj, dx, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JSparseB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) override {
        JSparseB_model_robertson(JB, t, x, p, k, h, cj, xB, dx, dxB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_Jrz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        Jrz_model_robertson(nllh, iz, p, k, rz, sigmaz);
    }

    virtual int model_Jv(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *v, const realtype *w, const realtype *dwdx) override {
        Jv_model_robertson(Jv, t, x, p, k, h, cj, dx, v, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual int model_JvB(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *vB, const realtype *w, const realtype *dwdx) override {
        JvB_model_robertson(JvB, t, x, p, k, h, cj, xB, dx, dxB, vB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_Jy(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        Jy_model_robertson(nllh, iy, p, k, y, sigmay, my);
    }

    virtual void model_Jz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        Jz_model_robertson(nllh, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_M(realtype *M, const realtype t, const realtype *x, const realtype *p, const realtype *k) override {
        M_model_robertson(M, t, x, p, k);
    }

    virtual void model_dJrzdsigma(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        dJrzdsigma_model_robertson(dJrzdsigma, iz, p, k, rz, sigmaz);
    }

    virtual void model_dJrzdz(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
        dJrzdz_model_robertson(dJrzdz, iz, p, k, rz, sigmaz);
    }

    virtual void model_dJydsigma(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydsigma_model_robertson(dJydsigma, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJydy(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydy_model_robertson(dJydy, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJzdsigma(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        dJzdsigma_model_robertson(dJzdsigma, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_dJzdz(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
        dJzdz_model_robertson(dJzdz, iz, p, k, z, sigmaz, mz);
    }

    virtual void model_deltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
        deltaqB_model_robertson(deltaqB, t, x, p, k, h, ip, ie, xdot, xdot_old, xB);
    }

    virtual void model_deltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) override {
        deltasx_model_robertson(deltasx, t, x, p, k, h, w, ip, ie, xdot, xdot_old, sx, stau);
    }

    virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {
        deltax_model_robertson(deltax, t, x, p, k, h, ie, xdot, xdot_old);
    }

    virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
        deltaxB_model_robertson(deltaxB, t, x, p, k, h, ie, xdot, xdot_old, xB);
    }

    virtual void model_drzdp(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        drzdp_model_robertson(drzdp, ie, t, x, p, k, h, ip);
    }

    virtual void model_drzdx(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        drzdx_model_robertson(drzdx, ie, t, x, p, k, h);
    }

    virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
        dsigma_ydp_model_robertson(dsigmaydp, t, p, k, ip);
    }

    virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
        dsigma_zdp_model_robertson(dsigmazdp, t, p, k, ip);
    }

    virtual void model_dwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdp_model_robertson(dwdp, t, x, p, k, h, w);
    }

    virtual void model_dwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdx_model_robertson(dwdx, t, x, p, k, h, w);
    }

    virtual int model_dxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *dx, const realtype *w, const realtype *dwdp) override {
        dxdotdp_model_robertson(dxdotdp, t, x, p, k, h, ip, dx, w, dwdp);
        return AMICI_SUCCESS;
    }

    virtual void model_dydp(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        dydp_model_robertson(dydp, t, x, p, k, h, ip);
    }

    virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dydx_model_robertson(dydx, t, x, p, k, h);
    }

    virtual void model_dzdp(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        dzdp_model_robertson(dzdp, ie, t, x, p, k, h, ip);
    }

    virtual void model_dzdx(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dzdx_model_robertson(dzdx, ie, t, x, p, k, h);
    }

    virtual int model_qBdot(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdp) override {
        qBdot_model_robertson(qBdot, ip, t, x, p, k, h, xB, dx, dxB, w, dwdp);
        return AMICI_SUCCESS;
    }

    virtual int model_root(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx) override {
        root_model_robertson(root, t, x, p, k, h, dx);
        return AMICI_SUCCESS;
    }

    virtual void model_rz(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        rz_model_robertson(rz, ie, t, x, p, k, h);
    }

    virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) override {
        sigma_y_model_robertson(sigmay, t, p, k);
    }

    virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) override {
        sigma_z_model_robertson(sigmaz, t, p, k);
    }

    virtual void model_srz(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
        srz_model_robertson(srz, ie, t, x, p, k, h, sx, ip);
    }

    virtual void model_stau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) override {
        stau_model_robertson(stau, t, x, p, k, h, sx, ip, ie);
    }

    virtual void model_sx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) override {
        sx0_model_robertson(sx0, t, x0, p, k, ip);
    }

    virtual int model_sxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *dx, const realtype *sx, const realtype *sdx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *M, const realtype *dxdotdp) override {
        sxdot_model_robertson(sxdot, t, x, p, k, h, ip, dx, sx, sdx, w, dwdx, J, M, dxdotdp);
        return AMICI_SUCCESS;
    }

    virtual void model_sz(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
        sz_model_robertson(sz, ie, t, x, p, k, h, sx, ip);
    }

    virtual void model_w(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        w_model_robertson(w, t, x, p, k, h);
    }

    virtual void model_x0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_model_robertson(x0, t, p, k);
    }

    virtual int model_xBdot(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) override {
        xBdot_model_robertson(xBdot, t, x, p, k, h, xB, dx, dxB, w, dwdx);
        return AMICI_SUCCESS;
    }

    virtual void model_xdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) override {
        xdot_model_robertson(xdot, t, x, p, k, h, dx, w);
    }

    virtual void model_y(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        y_model_robertson(y, t, x, p, k, h);
    }

    virtual void model_z(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        z_model_robertson(z, ie, t, x, p, k, h);
    }

};

#endif /* _amici_wrapfunctions_h */
