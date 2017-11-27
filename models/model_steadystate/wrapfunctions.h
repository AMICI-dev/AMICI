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
extern void Jv_model_steadystate(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx);
extern void JvB_model_steadystate(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx);
extern void Jy_model_steadystate(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydsigma_model_steadystate(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydy_model_steadystate(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dwdp_model_steadystate(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dwdx_model_steadystate(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_model_steadystate(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp);
extern void dydx_model_steadystate(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void qBdot_model_steadystate(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp);
extern void sigma_y_model_steadystate(double *sigmay, const realtype t, const realtype *p, const realtype *k);
extern void sxdot_model_steadystate(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp);
extern void w_model_steadystate(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void x0_model_steadystate(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void xBdot_model_steadystate(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void xdot_model_steadystate(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_model_steadystate(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

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
    }

    virtual void model_dJrzdsigma(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
    }

    virtual void model_dJrzdz(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
    }

    virtual void model_dJydsigma(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydsigma_model_steadystate(dJydsigma, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJydy(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydy_model_steadystate(dJydy, iy, p, k, y, sigmay, my);
    }

    virtual void model_dJzdsigma(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
    }

    virtual void model_dJzdz(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
    }

    virtual void model_deltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
    }

    virtual void model_deltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) override {
    }

    virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {
    }

    virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
    }

    virtual void model_drzdp(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
    }

    virtual void model_drzdx(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
    }

    virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
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
    }

    virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dydx_model_steadystate(dydx, t, x, p, k, h);
    }

    virtual void model_dzdp(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
    }

    virtual void model_dzdx(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    virtual int model_qBdot(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) override {
        qBdot_model_steadystate(qBdot, ip, t, x, p, k, h, xB, w, dwdp);
        return AMICI_SUCCESS;
    }

    virtual int model_root(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        return AMICI_SUCCESS;
    }

    virtual void model_rz(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) override {
        sigma_y_model_steadystate(sigmay, t, p, k);
    }

    virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) override {
    }

    virtual void model_srz(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
    }

    virtual void model_stau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) override {
    }

    virtual void model_sx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) override {
    }

    virtual int model_sxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) override {
        sxdot_model_steadystate(sxdot, t, x, p, k, h, ip, sx, w, dwdx, J, dxdotdp);
        return AMICI_SUCCESS;
    }

    virtual void model_sz(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
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
    }

};

#endif /* _amici_wrapfunctions_h */
