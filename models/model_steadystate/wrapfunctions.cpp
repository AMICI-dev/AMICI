#include "wrapfunctions.h"

#include <include/cvodewrap.h>

UserData getUserData(){
    return UserData(5,
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
                    AMICI_SCALING_LOG10,
                    AMICI_O2MODE_NONE);
}

int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t){
   return CVodeInit(cvode_mem, xdot_model_steadystate, RCONST(t), x);
}

int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t){
    return CVodeInitB(cvode_mem, which, xBdot_model_steadystate, RCONST(t), xB);
}

int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot){
    return CVodeQuadInitB(cvode_mem, which, qBdot_model_steadystate, qBdot);
}

int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data){
                    UserData *udata = (UserData*) user_data;
    return CVodeSensInit1(cvode_mem, udata->nplist, udata->sensi_meth, sxdot_model_steadystate, sx);
}

int wrap_RootInit(void *cvode_mem, void *user_data){
                    UserData *udata = (UserData*) user_data;
    return CVodeRootInit(cvode_mem, 0, root_model_steadystate);
}

int wrap_SetDenseJacFn(void *cvode_mem){
    return CVDlsSetDenseJacFn(cvode_mem, J_model_steadystate);
}

int wrap_SetSparseJacFn(void *cvode_mem){
    return CVSlsSetSparseJacFn(cvode_mem, JSparse_model_steadystate);
}

int wrap_SetBandJacFn(void *cvode_mem){
    return CVDlsSetBandJacFn(cvode_mem, JBand_model_steadystate);
}

int wrap_SetJacTimesVecFn(void *cvode_mem){
    return CVSpilsSetJacTimesVecFn(cvode_mem, Jv_model_steadystate);
}

int wrap_SetDenseJacFnB(void *cvode_mem,int which){
    return CVDlsSetDenseJacFnB(cvode_mem, which, JB_model_steadystate);
}

int wrap_SetSparseJacFnB(void *cvode_mem,int which){
    return CVSlsSetSparseJacFnB(cvode_mem, which, JSparseB_model_steadystate);
}

int wrap_SetBandJacFnB(void *cvode_mem,int which){
    return CVDlsSetBandJacFnB(cvode_mem, which, JBandB_model_steadystate);
}

int wrap_SetJacTimesVecFnB(void *cvode_mem,int which){
    return CVSpilsSetJacTimesVecFnB(cvode_mem, which, JvB_model_steadystate);
}

int fx0(N_Vector x0, void *user_data){
    return x0_model_steadystate(x0, user_data);
}

int fdx0(N_Vector x0, N_Vector dx0, void *user_data){
    return 0;
}

int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data){
    return sx0_model_steadystate(sx0, x, dx, user_data);
}

int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data){
    return 0;
}

int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    return J_model_steadystate(N, t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
    return JB_model_steadystate(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

int fJDiag(realtype t, N_Vector JDiag, N_Vector x, void *user_data){
    return JDiag_model_steadystate(t, JDiag, x, user_data);
}

int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp){
    return Jv_model_steadystate(v, Jv, t, x, xdot, user_data, tmp);
}

int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data){
    return root_model_steadystate(t, x, root, user_data);
}

int frz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata){
    return rz_model_steadystate(t, ie, x, user_data, tdata, rdata);
}

int fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata){
    return srz_model_steadystate(t, ie, x, sx, user_data, tdata, rdata);
}

int fstau(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata){
    return stau_model_steadystate(t, ie, x, sx, user_data, tdata);
}

int fy(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata){
    return y_model_steadystate(t, it, x, user_data, rdata);
}

int fdydp(realtype t, int it, N_Vector x, void *user_data, TempData *tdata){
    return dydp_model_steadystate(t, it, x, user_data, tdata);
}

int fdydx(realtype t, int it, N_Vector x, void *user_data, TempData *tdata){
    return dydx_model_steadystate(t, it, x, user_data, tdata);
}

int fz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata){
    return z_model_steadystate(t, ie, x, user_data, tdata, rdata);
}

int fsz(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata){
    return sz_model_steadystate(t, ie, x, sx, user_data, tdata, rdata);
}

int fdzdp(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata){
    return dzdp_model_steadystate(t, ie, x, user_data, tdata);
}

int fdzdx(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata){
    return dzdx_model_steadystate(t, ie, x, user_data, tdata);
}

int fdrzdp(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata){
    return drzdp_model_steadystate(t, ie, x, user_data, tdata);
}

int fdrzdx(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata){
    return drzdx_model_steadystate(t, ie, x, user_data, tdata);
}

int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data){
    return xdot_model_steadystate(t, x, xdot, user_data);
}

int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data){
    return xBdot_model_steadystate(t, x, xB, xBdot, user_data);
}

int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data){
    return qBdot_model_steadystate(t, x, xB, qBdot, user_data);
}

int fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data){
    return dxdotdp_model_steadystate(t, x, dx, user_data);
}

int fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata){
    return deltax_model_steadystate(t, ie, x, xdot, xdot_old, user_data, tdata);
}

int fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data, TempData *tdata){
    return deltasx_model_steadystate(t, ie, x, xdot, xdot_old, sx, user_data, tdata);
}

int fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata){
    return deltaxB_model_steadystate(t, ie, x, xB, xdot, xdot_old, user_data, tdata);
}

int fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata){
    return deltaqB_model_steadystate(t, ie, x, xB, qBdot, xdot, xdot_old, user_data, tdata);
}

int fsigma_y(realtype t, void *user_data, TempData *tdata){
    return sigma_y_model_steadystate(t, user_data, tdata);
}

int fdsigma_ydp(realtype t, void *user_data, TempData *tdata){
    return dsigma_ydp_model_steadystate(t, user_data, tdata);
}

int fsigma_z(realtype t, int ie, void *user_data, TempData *tdata){
    return sigma_z_model_steadystate(t, ie, user_data, tdata);
}

int fdsigma_zdp(realtype t, int ie, void *user_data, TempData *tdata){
    return dsigma_zdp_model_steadystate(t, ie, user_data, tdata);
}

int fJy(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jy_model_steadystate(t, it, x, user_data, tdata, edata, rdata);
}

int fJz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jz_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

int fJrz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return Jrz_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

int fdJydy(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJydy_model_steadystate(t, it, x, user_data, tdata, edata, rdata);
}

int fdJydsigma(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJydsigma_model_steadystate(t, it, x, user_data, tdata, edata, rdata);
}

int fdJzdz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJzdz_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

int fdJzdsigma(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJzdsigma_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

int fdJrzdz(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJrzdz_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

int fdJrzdsigma(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata){
    return dJrzdsigma_model_steadystate(t, ie, x, user_data, tdata, edata, rdata);
}

