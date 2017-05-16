                
#include "wrapfunctions.h"
#include <include/udata_accessors.h>
                
                void init_modeldims(UserData *udata){
                   np = 17;
                   nx = 9;
                   nxtrue = 9;
                   nk = 2;
                   ny = 3;
                   nytrue = 3;
                   nz = 0;
                   nztrue = 0;
                   ne = 0;
                   ng = 1;
                   nw = 2;
                   ndwdx = 1;
                   ndwdp = 5;
                   nnz = 18;
                   ubw = 8;
                   lbw = 1;
                   udata->am_pscale = AMI_SCALING_LOG10;
                   udata->am_o2mode = AMI_O2MODE_NONE;
                }
                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t){
                    return CVodeInit(cvode_mem, xdot_model_jakstat_adjoint, RCONST(t), x);
                }
                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t){
                    return CVodeInitB(cvode_mem, which, xBdot_model_jakstat_adjoint, RCONST(t), xB);
                }
                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot){
                    return CVodeQuadInitB(cvode_mem, which, qBdot_model_jakstat_adjoint, qBdot);
                }
                int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data){
                    UserData *udata = (UserData*) user_data;
                    return CVodeSensInit1(cvode_mem, nplist, sensi_meth, sxdot_model_jakstat_adjoint, sx);
                }
                
                int wrap_RootInit(void *cvode_mem, void *user_data){
                    UserData *udata = (UserData*) user_data;
                    return CVodeRootInit(cvode_mem, 0, root_model_jakstat_adjoint);
                }
                
                int wrap_SetDenseJacFn(void *cvode_mem){
                    return CVDlsSetDenseJacFn(cvode_mem, J_model_jakstat_adjoint);
                }
                int wrap_SetSparseJacFn(void *cvode_mem){
                    return CVSlsSetSparseJacFn(cvode_mem, JSparse_model_jakstat_adjoint);
                }
                int wrap_SetBandJacFn(void *cvode_mem){
                    return CVDlsSetBandJacFn(cvode_mem, JBand_model_jakstat_adjoint);
                }
                int wrap_SetJacTimesVecFn(void *cvode_mem){
                    return CVSpilsSetJacTimesVecFn(cvode_mem, Jv_model_jakstat_adjoint);
                }
                int wrap_SetDenseJacFnB(void *cvode_mem,int which){
                    return CVDlsSetDenseJacFnB(cvode_mem, which, JB_model_jakstat_adjoint);
                }
                int wrap_SetSparseJacFnB(void *cvode_mem,int which){
                    return CVSlsSetSparseJacFnB(cvode_mem, which, JSparseB_model_jakstat_adjoint);
                }
                int wrap_SetBandJacFnB(void *cvode_mem,int which){
                    return CVDlsSetBandJacFnB(cvode_mem, which, JBandB_model_jakstat_adjoint);
                }
                int wrap_SetJacTimesVecFnB(void *cvode_mem,int which){
                    return CVSpilsSetJacTimesVecFnB(cvode_mem, which, JvB_model_jakstat_adjoint);
                }
                int fx0(N_Vector x0, void *user_data){
                    return x0_model_jakstat_adjoint(x0, user_data);
                }
                
                int fdx0(N_Vector x0, N_Vector dx0, void *user_data){
                    return(0);
                }
                
                int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data){
                    return sx0_model_jakstat_adjoint(sx0, x, dx, user_data);
                }
                
                int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data){
                    return(0);
                }
                
                int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
                    return J_model_jakstat_adjoint(N, t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
                }
                
                int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B){
                    return JB_model_jakstat_adjoint(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
                }
                
                int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data){
                    return root_model_jakstat_adjoint(t, x, root, user_data);
                }
                
                int fsroot(realtype t, int ie, int *nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data){
                    return sroot_model_jakstat_adjoint(t, ie, nroots, sroot, x, sx, user_data);
                }
                
                int fs2root(realtype t, int ie, int *nroots, realtype *s2root, N_Vector x, N_Vector *sx, void *user_data){
                    warnMsgIdAndTxt("AMICI:mex:s2root:NotAvailable","ERROR: The function s2root was called but not compiled for this model.");
                    return(-1);
                }
                
                int fstau(realtype t, int ie, realtype *stau, N_Vector x, N_Vector *sx, void *user_data){
                    return stau_model_jakstat_adjoint(t, ie, stau, x, sx, user_data);
                }
                
                int fy(realtype t, int it, realtype *y, N_Vector x, void *user_data){
                    return y_model_jakstat_adjoint(t, it, y, x, user_data);
                }
                
                int fdydp(realtype t, int it, realtype *dydp, N_Vector x, void *user_data){
                    return dydp_model_jakstat_adjoint(t, it, dydp, x, user_data);
                }
                
                int fdydx(realtype t, int it, realtype *dydx, N_Vector x, void *user_data){
                    return dydx_model_jakstat_adjoint(t, it, dydx, x, user_data);
                }
                
                int fz(realtype t, int ie, int *nroots, realtype *z, N_Vector x, void *user_data){
                    return z_model_jakstat_adjoint(t, ie, nroots, z, x, user_data);
                }
                
                int fsz(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data){
                    return sz_model_jakstat_adjoint(t, ie, nroots, sz, x, sx, user_data);
                }
                
                int fsz_tf(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data){
                    return sz_tf_model_jakstat_adjoint(t, ie, nroots, sz, x, sx, user_data);
                }
                
                int fdzdp(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data){
                    return dzdp_model_jakstat_adjoint(t, ie, dzdp, x, user_data);
                }
                
                int fdzdx(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data){
                    return dzdx_model_jakstat_adjoint(t, ie, dzdx, x, user_data);
                }
                
                int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data){
                    return xdot_model_jakstat_adjoint(t, x, xdot, user_data);
                }
                
                int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data){
                    return xBdot_model_jakstat_adjoint(t, x, xB, xBdot, user_data);
                }
                
                int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data){
                    return qBdot_model_jakstat_adjoint(t, x, xB, qBdot, user_data);
                }
                
                int fdxdotdp(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data){
                    return dxdotdp_model_jakstat_adjoint(t, dxdotdp, x, dx, user_data);
                }
                
                int fdeltax(realtype t, int ie, realtype *deltax, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data){
                    return deltax_model_jakstat_adjoint(t, ie, deltax, x, xdot, xdot_old, user_data);
                }
                
                int fdeltasx(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data){
                    return deltasx_model_jakstat_adjoint(t, ie, deltasx, x, xdot, xdot_old, sx, user_data);
                }
                
                int fdeltaxB(realtype t, int ie, realtype *deltaxB, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data){
                    return deltaxB_model_jakstat_adjoint(t, ie, deltaxB, x, xB, xdot, xdot_old, user_data);
                }
                
                int fdeltaqB(realtype t, int ie, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data){
                    return deltaqB_model_jakstat_adjoint(t, ie, deltaqB, x, xB, qBdot, xdot, xdot_old, user_data);
                }
                
                int fsigma_y(realtype t, realtype *sigma_y, void *user_data){
                    return sigma_y_model_jakstat_adjoint(t, sigma_y, user_data);
                }
                
                int fdsigma_ydp(realtype t, realtype *dsigma_ydp, void *user_data){
                    return dsigma_ydp_model_jakstat_adjoint(t, dsigma_ydp, user_data);
                }
                
                int fsigma_z(realtype t, int ie, realtype *sigma_z, void *user_data){
                    return sigma_z_model_jakstat_adjoint(t, ie, sigma_z, user_data);
                }
                
                int fdsigma_zdp(realtype t, int ie, realtype *dsigma_zdp, void *user_data){
                    return dsigma_zdp_model_jakstat_adjoint(t, ie, dsigma_zdp, user_data);
                }
                
                int fJy(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data){
                    return Jy_model_jakstat_adjoint(t, it, Jy, y, x, my, sigma_y, user_data);
                }
                
                int fJz(realtype t, int ie, realtype *Jz, realtype *z, N_Vector x, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data){
                    return Jz_model_jakstat_adjoint(t, ie, Jz, z, x, mz, sigma_z, user_data, temp_data);
                }
                
                int fdJydx(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data){
                    return dJydx_model_jakstat_adjoint(t, it, dJydx, y, x, dydx, my, sigma_y, user_data);
                }
                
                int fdJydy(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data){
                    return dJydy_model_jakstat_adjoint(t, it, dJydy, y, x, my, sigma_y, user_data);
                }
                
                int fdJydp(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data){
                    return dJydp_model_jakstat_adjoint(t, it, dJydp, y, x, dydp, my, sigma_y, dsigma_ydp, user_data);
                }
                
                int fdJzdx(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data){
                    return dJzdx_model_jakstat_adjoint(t, ie, dJzdx, z, x, dzdx, mz, sigma_z, user_data, temp_data);
                }
                
                int fdJzdp(realtype t, int ie, realtype *dJzdp, realtype *z, N_Vector x, realtype *dzdp, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data){
                    return dJzdp_model_jakstat_adjoint(t, ie, dJzdp, z, x, dzdp, mz, sigma_z, dsigma_zdp, user_data, temp_data);
                }
                
                int fsJz(realtype t, int ie, realtype *sJz, realtype *s2Jz, realtype *dJzdz, realtype *dJzdp, realtype *sz, realtype *dzdp, realtype *mz, void *user_data, void *temp_data){
                    return sJz_model_jakstat_adjoint(t, ie, sJz, s2Jz, dJzdz, dJzdp, sz, dzdp, mz, user_data, temp_data);
                }
                
