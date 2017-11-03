#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h
#include <math.h>
#include <include/amici_model.h>
#include "model_nested_events.h"
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
void fx0(N_Vector x0, void *user_data);
void fdx0(N_Vector x0, N_Vector dx0, void *user_data);
void fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data);
void fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data);
void fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data);
void fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
void froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data);
void frz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata);
void fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata);
void fstau(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata);
void fy(realtype t, int it, N_Vector x, void *user_data, amici::ReturnData *rdata);
void fdydp(realtype t, int it, N_Vector x, amici::TempData *tdata);
void fdydx(realtype t, int it, N_Vector x, amici::TempData *tdata);
void fz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata);
void fsz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata);
void fdzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata);
void fdzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata);
void fdrzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata);
void fdrzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata);
void fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data);
void fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data);
void fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data);
void fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data);
void fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata);
void fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, amici::TempData *tdata);
void fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata);
void fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata);
void fsigma_y(realtype t, amici::TempData *tdata);
void fdsigma_ydp(realtype t, amici::TempData *tdata);
void fsigma_z(realtype t, int ie, amici::TempData *tdata);
void fdsigma_zdp(realtype t, int ie, amici::TempData *tdata);
void fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
void fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
void fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2);
void fJy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fJz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fJrz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJydy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJydsigma(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJrzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);
void fdJrzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata);

class Model_model_nested_events : public amici::Model {
public:
    Model_model_nested_events() : amici::Model(5,
                    1,
                    1,
                    0,
                    1,
                    1,
                    0,
                    0,
                    3,
                    1,
                    0,
                    0,
                    0,
                    1,
                    0,
                    0,
                    amici::AMICI_O2MODE_NONE)
    {
        z2event = new int[0] {};
        idlist = new realtype[1] {0,};
    }

    amici::Solver *getSolver() override {
        return new amici::CVodeSolver();
    }

    int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) override {
        J_model_nested_events(N, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
        return(0);
    }

    int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) override {
        JB_model_nested_events(NeqBdot, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
        return(0);
    }

    int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) override {
        JBand_model_nested_events(N, mupper, mlower, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
        return(0);
    }

    int fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) override {
        JBandB_model_nested_events(NeqBdot, mupper, mlower, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
        return(0);
    }

    void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data) override {
        JDiag_model_nested_events(t, JDiag, cj, x, dx, user_data);
        return;
    }

    int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) override {
        JSparse_model_nested_events(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
        return(0);
    }

    int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) override {
        JSparseB_model_nested_events(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
        return(0);
    }

    void fJrz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        Jrz_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) override {
        Jv_model_nested_events(t, x, dx, xdot, v, Jv, cj, user_data, tmp1, tmp2);
        return(0);
    }

    int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2) override {
        JvB_model_nested_events(t, x, dx, xB, dxB, xBdot, vB, JvB, cj, user_data, tmpB1, tmpB2);
        return(0);
    }

    void fJy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        Jy_model_nested_events(t, it, x, tdata, edata, rdata);
        return;
    }

    void fJz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        Jz_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    void fdJrzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJrzdsigma_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    void fdJrzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJrzdz_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    void fdJydsigma(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJydsigma_model_nested_events(t, it, x, tdata, edata, rdata);
        return;
    }

    void fdJydy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJydy_model_nested_events(t, it, x, tdata, edata, rdata);
        return;
    }

    void fdJzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJzdsigma_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    void fdJzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) override {
        dJzdz_model_nested_events(t, ie, x, tdata, edata, rdata);
        return;
    }

    void fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata) override {
        deltaqB_model_nested_events(t, ie, x, xB, qBdot, xdot, xdot_old, tdata);
        return;
    }

    void fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, amici::TempData *tdata) override {
        deltasx_model_nested_events(t, ie, x, xdot, xdot_old, sx, tdata);
        return;
    }

    void fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata) override {
        deltax_model_nested_events(t, ie, x, xdot, xdot_old, tdata);
        return;
    }

    void fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata) override {
        deltaxB_model_nested_events(t, ie, x, xB, xdot, xdot_old, tdata);
        return;
    }

    void fdrzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata) override {
        drzdp_model_nested_events(t, ie, x, tdata);
        return;
    }

    void fdrzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata) override {
        drzdx_model_nested_events(t, ie, x, tdata);
        return;
    }

    void fdsigma_ydp(realtype t, amici::TempData *tdata) override {
        dsigma_ydp_model_nested_events(t, tdata);
        return;
    }

    void fdsigma_zdp(realtype t, int ie, amici::TempData *tdata) override {
        dsigma_zdp_model_nested_events(t, ie, tdata);
        return;
    }

    void fdwdp(realtype t, N_Vector x, N_Vector dx, void *user_data) override {
        dwdp_model_nested_events(t, x, dx, user_data);
        return;
    }

    void fdwdx(realtype t, N_Vector x, N_Vector dx, void *user_data) override {
        dwdx_model_nested_events(t, x, dx, user_data);
        return;
    }

    void fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data) override {
        dxdotdp_model_nested_events(t, x, dx, user_data);
        return;
    }

    void fdydp(realtype t, int it, N_Vector x, amici::TempData *tdata) override {
        dydp_model_nested_events(t, it, x, tdata);
        return;
    }

    void fdydx(realtype t, int it, N_Vector x, amici::TempData *tdata) override {
        dydx_model_nested_events(t, it, x, tdata);
        return;
    }

    void fdzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata) override {
        dzdp_model_nested_events(t, ie, x, tdata);
        return;
    }

    void fdzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata) override {
        dzdx_model_nested_events(t, ie, x, tdata);
        return;
    }

    int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data) override {
        qBdot_model_nested_events(t, x, dx, xB, dxB, qBdot, user_data);
        return(0);
    }

    int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) override {
        root_model_nested_events(t, x, dx, root, user_data);
        return(0);
    }

    void frz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata) override {
        rz_model_nested_events(t, ie, x, tdata, rdata);
        return;
    }

    void fsigma_y(realtype t, amici::TempData *tdata) override {
        sigma_y_model_nested_events(t, tdata);
        return;
    }

    void fsigma_z(realtype t, int ie, amici::TempData *tdata) override {
        sigma_z_model_nested_events(t, ie, tdata);
        return;
    }

    void fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata) override {
        srz_model_nested_events(t, ie, x, sx, tdata, rdata);
        return;
    }

    void fstau(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata) override {
        stau_model_nested_events(t, ie, x, sx, tdata);
        return;
    }

    void fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) override {
        sx0_model_nested_events(sx0, x, dx, user_data);
        return;
    }

    int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) override {
        sxdot_model_nested_events(Ns, t, x, dx, xdot, ip, sx, sdx, sxdot, user_data, tmp1, tmp2, tmp3);
        return(0);
    }

    void fsz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata) override {
        sz_model_nested_events(t, ie, x, sx, tdata, rdata);
        return;
    }

    void fw(realtype t, N_Vector x, N_Vector dx, void *user_data) override {
        w_model_nested_events(t, x, dx, user_data);
        return;
    }

    void fx0(N_Vector x0, void *user_data) override {
        x0_model_nested_events(x0, user_data);
        return;
    }

    int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data) override {
        xBdot_model_nested_events(t, x, dx, xB, dxB, xBdot, user_data);
        return(0);
    }

    int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data) override {
        xdot_model_nested_events(t, x, dx, xdot, user_data);
        return(0);
    }

    void fy(realtype t, int it, N_Vector x, void *user_data, amici::ReturnData *rdata) override {
        y_model_nested_events(t, it, x, user_data, rdata);
        return;
    }

    void fz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata) override {
        z_model_nested_events(t, ie, x, tdata, rdata);
        return;
    }

};

#endif /* _amici_wrapfunctions_h */
