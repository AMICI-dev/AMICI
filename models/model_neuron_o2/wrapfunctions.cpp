#include "wrapfunctions.h"
#include <include/amici_model.h>
#include <include/udata.h>

Model *getModel() {
    return new Model_model_neuron_o2();
}

int fx0(N_Vector x0, void *user_data){
    return x0_model_neuron_o2(x0, user_data);
}

int fdx0(N_Vector x0, N_Vector dx0, void *user_data){
    return 0;
}

int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data){
    return sx0_model_neuron_o2(sx0, x, dx, user_data);
}

int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data){
    return 0;
}

int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return J_model_neuron_o2(N, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JB_model_neuron_o2(NeqBdot, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

int fJDiag(realtype t, N_Vector JDiag, realtype cj, N_Vector x, N_Vector dx, void *user_data){
    return JDiag_model_neuron_o2(t, JDiag, cj, x, dx, user_data);
}

int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2){
    return Jv_model_neuron_o2(t, x, dx, xdot, v, Jv, cj, user_data, tmp1, tmp2);
}

int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data){
    return root_model_neuron_o2(t, x, dx, root, user_data);
}

int frz(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata){
    return rz_model_neuron_o2(t, ie, x, tdata, rdata);
}

int fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata){
    return srz_model_neuron_o2(t, ie, x, sx, tdata, rdata);
}

int fstau(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata){
    return stau_model_neuron_o2(t, ie, x, sx, tdata);
}

int fy(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata){
    return y_model_neuron_o2(t, it, x, user_data, rdata);
}

int fdydp(realtype t, int it, N_Vector x, TempData *tdata){
    return dydp_model_neuron_o2(t, it, x, tdata);
}

int fdydx(realtype t, int it, N_Vector x, TempData *tdata){
    return dydx_model_neuron_o2(t, it, x, tdata);
}

int fz(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata){
    return z_model_neuron_o2(t, ie, x, tdata, rdata);
}

int fsz(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata){
    return sz_model_neuron_o2(t, ie, x, sx, tdata, rdata);
}

int fdzdp(realtype t, int ie, N_Vector x, TempData *tdata){
    return dzdp_model_neuron_o2(t, ie, x, tdata);
}

int fdzdx(realtype t, int ie, N_Vector x, TempData *tdata){
    return dzdx_model_neuron_o2(t, ie, x, tdata);
}

int fdrzdp(realtype t, int ie, N_Vector x, TempData *tdata){
    return drzdp_model_neuron_o2(t, ie, x, tdata);
}

int fdrzdx(realtype t, int ie, N_Vector x, TempData *tdata){
    return drzdx_model_neuron_o2(t, ie, x, tdata);
}

int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,int ip,  N_Vector sx, N_Vector sdx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return sxdot_model_neuron_o2(Ns, t, x, dx, xdot, ip, sx, sdx, sxdot, user_data, tmp1, tmp2, tmp3);
}

int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data){
    return xdot_model_neuron_o2(t, x, dx, xdot, user_data);
}

int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, void *user_data){
    return xBdot_model_neuron_o2(t, x, dx, xB, dxB, xBdot, user_data);
}

int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data){
    return qBdot_model_neuron_o2(t, x, dx, xB, dxB, qBdot, user_data);
}

int fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data){
    return dxdotdp_model_neuron_o2(t, x, dx, user_data);
}

int fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, TempData *tdata){
    return deltax_model_neuron_o2(t, ie, x, xdot, xdot_old, tdata);
}

int fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, TempData *tdata){
    return deltasx_model_neuron_o2(t, ie, x, xdot, xdot_old, sx, tdata);
}

int fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, TempData *tdata){
    return deltaxB_model_neuron_o2(t, ie, x, xB, xdot, xdot_old, tdata);
}

int fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, TempData *tdata){
    return deltaqB_model_neuron_o2(t, ie, x, xB, qBdot, xdot, xdot_old, tdata);
}

int fsigma_y(realtype t, TempData *tdata){
    return sigma_y_model_neuron_o2(t, tdata);
}

int fdsigma_ydp(realtype t, TempData *tdata){
    return dsigma_ydp_model_neuron_o2(t, tdata);
}

int fsigma_z(realtype t, int ie, TempData *tdata){
    return sigma_z_model_neuron_o2(t, ie, tdata);
}

int fdsigma_zdp(realtype t, int ie, TempData *tdata){
    return dsigma_zdp_model_neuron_o2(t, ie, tdata);
}

int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return JSparse_model_neuron_o2(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return JBand_model_neuron_o2(N, mupper, mlower, t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JSparseB_model_neuron_o2(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

int fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JBandB_model_neuron_o2(NeqBdot, mupper, mlower, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2){
    return JvB_model_neuron_o2(t, x, dx, xB, dxB, xBdot, vB, JvB, cj, user_data, tmpB1, tmpB2);
}

int fJy(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jy_model_neuron_o2(t, it, x, tdata, edata, rdata);
}

int fJz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jz_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

int fJrz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jrz_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

int fdJydy(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJydy_model_neuron_o2(t, it, x, tdata, edata, rdata);
}

int fdJydsigma(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJydsigma_model_neuron_o2(t, it, x, tdata, edata, rdata);
}

int fdJzdz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJzdz_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

int fdJzdsigma(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJzdsigma_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

int fdJrzdz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJrzdz_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

int fdJrzdsigma(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJrzdsigma_model_neuron_o2(t, ie, x, tdata, edata, rdata);
}

