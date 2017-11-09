#include "wrapfunctions.h"
#include <include/amici_model.h>
#include <include/udata.h>

using namespace amici;

Model *getModel() {
    return new Model_model_events();
}

void fx0(N_Vector x0, void *user_data){
    return x0_model_events(x0, user_data);
}

void fdx0(N_Vector x0, N_Vector dx0, void *user_data){
    return;
}

void fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data){
    return sx0_model_events(sx0, x, dx, user_data);
}

void fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data){
    return;
}

void fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return J_model_events(N, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

void fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JB_model_events(NeqBdot, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

void fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data){
    return JDiag_model_events(t, JDiag, cj, x, dx, user_data);
}

void fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2){
    return Jv_model_events(t, x, dx, xdot, v, Jv, cj, user_data, tmp1, tmp2);
}

void froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data){
    return root_model_events(t, x, dx, root, user_data);
}

void frz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata){
    return rz_model_events(t, ie, x, tdata, rdata);
}

void fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata){
    return srz_model_events(t, ie, x, sx, tdata, rdata);
}

void fstau(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata){
    return stau_model_events(t, ie, x, sx, tdata);
}

void fy(realtype t, int it, N_Vector x, void *user_data, amici::ReturnData *rdata){
    return y_model_events(t, it, x, user_data, rdata);
}

void fdydp(realtype t, int it, N_Vector x, amici::TempData *tdata){
    return dydp_model_events(t, it, x, tdata);
}

void fdydx(realtype t, int it, N_Vector x, amici::TempData *tdata){
    return dydx_model_events(t, it, x, tdata);
}

void fz(realtype t, int ie, N_Vector x, amici::TempData *tdata, amici::ReturnData *rdata){
    return z_model_events(t, ie, x, tdata, rdata);
}

void fsz(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata){
    return sz_model_events(t, ie, x, sx, tdata, rdata);
}

void fdzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata){
    return dzdp_model_events(t, ie, x, tdata);
}

void fdzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata){
    return dzdx_model_events(t, ie, x, tdata);
}

void fdrzdp(realtype t, int ie, N_Vector x, amici::TempData *tdata){
    return drzdp_model_events(t, ie, x, tdata);
}

void fdrzdx(realtype t, int ie, N_Vector x, amici::TempData *tdata){
    return drzdx_model_events(t, ie, x, tdata);
}

void fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return sxdot_model_events(Ns, t, x, dx, xdot, ip, sx, sdx, sxdot, user_data, tmp1, tmp2, tmp3);
}

void fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data){
    return xdot_model_events(t, x, dx, xdot, user_data);
}

void fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data){
    return xBdot_model_events(t, x, dx, xB, dxB, xBdot, user_data);
}

void fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data){
    return qBdot_model_events(t, x, dx, xB, dxB, qBdot, user_data);
}

void fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data){
    return dxdotdp_model_events(t, x, dx, user_data);
}

void fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata){
    return deltax_model_events(t, ie, x, xdot, xdot_old, tdata);
}

void fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, amici::TempData *tdata){
    return deltasx_model_events(t, ie, x, xdot, xdot_old, sx, tdata);
}

void fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata){
    return deltaxB_model_events(t, ie, x, xB, xdot, xdot_old, tdata);
}

void fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, amici::TempData *tdata){
    return deltaqB_model_events(t, ie, x, xB, qBdot, xdot, xdot_old, tdata);
}

void fsigma_y(realtype t, amici::TempData *tdata){
    return sigma_y_model_events(t, tdata);
}

void fdsigma_ydp(realtype t, amici::TempData *tdata){
    return dsigma_ydp_model_events(t, tdata);
}

void fsigma_z(realtype t, int ie, amici::TempData *tdata){
    return sigma_z_model_events(t, ie, tdata);
}

void fdsigma_zdp(realtype t, int ie, amici::TempData *tdata){
    return dsigma_zdp_model_events(t, ie, tdata);
}

void fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return JSparse_model_events(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

void fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return JBand_model_events(N, mupper, mlower, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

void fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JSparseB_model_events(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

void fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JBandB_model_events(NeqBdot, mupper, mlower, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

void fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2){
    return JvB_model_events(t, x, dx, xB, dxB, xBdot, vB, JvB, cj, user_data, tmpB1, tmpB2);
}

void fJy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return Jy_model_events(t, it, x, tdata, edata, rdata);
}

void fJz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return Jz_model_events(t, ie, x, tdata, edata, rdata);
}

void fJrz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return Jrz_model_events(t, ie, x, tdata, edata, rdata);
}

void fdJydy(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJydy_model_events(t, it, x, tdata, edata, rdata);
}

void fdJydsigma(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJydsigma_model_events(t, it, x, tdata, edata, rdata);
}

void fdJzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJzdz_model_events(t, ie, x, tdata, edata, rdata);
}

void fdJzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJzdsigma_model_events(t, ie, x, tdata, edata, rdata);
}

void fdJrzdz(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJrzdz_model_events(t, ie, x, tdata, edata, rdata);
}

void fdJrzdsigma(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata){
    return dJrzdsigma_model_events(t, ie, x, tdata, edata, rdata);
}

