#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include <math.h>
#include <include/amici_model.h>
#include "model_events.h"
#include <include/amici_solver_cvodes.h>

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

    void J_model_events(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    model_J = &J_model_events;

    void JB_model_events(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    model_JB = &JB_model_events;

    void JBand_model_events(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    model_JBand = &JBand_model_events;

    void JBandB_model_events(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    model_JBandB = &JBandB_model_events;

    void JDiag_model_events(realtype t, N_Vector JDiag, N_Vector x, void *user_data);
    model_JDiag = &JDiag_model_events;

    void JSparse_model_events(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    model_JSparse = &JSparse_model_events;

    void JSparseB_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    model_JSparseB = &JSparseB_model_events;

    void Jrz_model_events(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz);
    model_Jrz = &Jrz_model_events;

    void Jv_model_events(realtype t, N_Vector x, N_Vector xdot, N_Vector v, N_Vector Jv, void *user_data, N_Vector tmp1, N_Vector tmp2);
    model_Jv = &Jv_model_events;

    void JvB_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, N_Vector vB, N_Vector JvB, void *user_data, N_Vector tmpB1, N_Vector tmpB2);
    model_JvB = &JvB_model_events;

    void Jy_model_events(double *nllh, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    model_Jy = &Jy_model_events;

    void Jz_model_events(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    model_Jz = &Jz_model_events;

    void dJrzdsigma_model_events(double *dJrzdsigma, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
    model_dJrzdsigma = &dJrzdsigma_model_events;

    void dJrzdz_model_events(double *dJrzdz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz);
    model_dJrzdz = &dJrzdz_model_events;

    void dJydsigma_model_events(double *dJydsigma, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    model_dJydsigma = &dJydsigma_model_events;

    void dJydy_model_events(double *dJydy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
    model_dJydy = &dJydy_model_events;

    void dJzdsigma_model_events(double *dJzdsigma, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    model_dJzdsigma = &dJzdsigma_model_events;

    void dJzdz_model_events(double *dJzdz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz);
    model_dJzdz = &dJzdz_model_events;

    void deltaqB_model_events(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB, const realtype *qBdot);
    model_deltaqB = &deltaqB_model_events;

    void deltasx_model_events(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau);
    model_deltasx = &deltasx_model_events;

    void deltax_model_events(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ie, const realtype *xdot, const realtype *xdot_old);
    model_deltax = &deltax_model_events;

    void deltaxB_model_events(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB);
    model_deltaxB = &deltaxB_model_events;

    void drzdp_model_events(double *drzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
    model_drzdp = &drzdp_model_events;

    void drzdx_model_events(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_drzdx = &drzdx_model_events;

    void dsigma_ydp_model_events(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip);
    model_dsigma_ydp = &dsigma_ydp_model_events;

    void dsigma_zdp_model_events(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip);
    model_dsigma_zdp = &dsigma_zdp_model_events;

    void dwdp_model_events(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *w);
    model_dwdp = &dwdp_model_events;

    void dwdx_model_events(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *w);
    model_dwdx = &dwdx_model_events;

    void dxdotdp_model_events(realtype t, N_Vector x,, void *user_data);
    model_dxdotdp = &dxdotdp_model_events;

    void dydp_model_events(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
    model_dydp = &dydp_model_events;

    void dydx_model_events(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_dydx = &dydx_model_events;

    void dzdp_model_events(double *dzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip);
    model_dzdp = &dzdp_model_events;

    void dzdx_model_events(double *dzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_dzdx = &dzdx_model_events;

    void qBdot_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data);
    model_qBdot = &qBdot_model_events;

    void root_model_events(realtype t, N_Vector x, realtype *root, void *user_data);
    model_root = &root_model_events;

    void rz_model_events(double *rz, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_rz = &rz_model_events;

    void sigma_y_model_events(double *sigmay, const realtype t, const realtype *p, const realtype *k);
    model_sigma_y = &sigma_y_model_events;

    void sigma_z_model_events(double *sigmaz, const realtype t, const realtype *p, const realtype *k);
    model_sigma_z = &sigma_z_model_events;

    void srz_model_events(double *srz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip);
    model_srz = &srz_model_events;

    void stau_model_events(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip, const int ie);
    model_stau = &stau_model_events;

    void sx0_model_events(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip);
    model_sx0 = &sx0_model_events;

    void sxdot_model_events(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2);
    model_sxdot = &sxdot_model_events;

    void sz_model_events(double *sz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip);
    model_sz = &sz_model_events;

    void w_model_events(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_w = &w_model_events;

    void x0_model_events(realtype *x0, const realtype t, const realtype *p, const realtype *k);
    model_x0 = &x0_model_events;

    void xBdot_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data);
    model_xBdot = &xBdot_model_events;

    void xdot_model_events(realtype t, N_Vector x, N_Vector xdot, void *user_data);
    model_xdot = &xdot_model_events;

    void y_model_events(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_y = &y_model_events;

    void z_model_events(double *z, const realtype t, const realtype *x, const realtype *p, const realtype *k);
    model_z = &z_model_events;

    amici::Solver *getSolver() override {
        return new amici::CVodeSolver();
    }

#endif /* _amici_wrapfunctions_h */
