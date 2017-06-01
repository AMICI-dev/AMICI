#ifndef _am_wrapfunctions_h
#define _am_wrapfunctions_h
#include <math.h>
#ifndef AMICI_WITHOUT_MATLAB
#include <mex.h>
#endif

#include <include/cvodewrap.h>

#include "model_jakstat_adjoint.h"

#include <include/udata.h>


#define pi M_PI

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC         UserData *getUserData();
                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t);
                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t);
                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot);
                int wrap_RootInit(void *cvode_mem, void *user_data);
                int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data);
                int wrap_SetDenseJacFn(void *cvode_mem);
                int wrap_SetSparseJacFn(void *cvode_mem);
                int wrap_SetBandJacFn(void *cvode_mem);
                int wrap_SetJacTimesVecFn(void *cvode_mem);
                int wrap_SetDenseJacFnB(void *cvode_mem,int which);
                int wrap_SetSparseJacFnB(void *cvode_mem,int which);
                int wrap_SetBandJacFnB(void *cvode_mem,int which);
                int wrap_SetJacTimesVecFnB(void *cvode_mem,int which);
                int fx0(N_Vector x0, void *user_data);
                int fdx0(N_Vector x0, N_Vector dx0, void *user_data);
                int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data);
                int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data);
                int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
                int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
                int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data);
                int fsroot(realtype t, int ie, int *nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data);
                int fs2root(realtype t, int ie, int *nroots, realtype *s2root, N_Vector x, N_Vector *sx, void *user_data);
                int fstau(realtype t, int ie, realtype *stau, N_Vector x, N_Vector *sx, void *user_data);
                int fy(realtype t, int it, realtype *y, N_Vector x, void *user_data);
                int fdydp(realtype t, int it, realtype *dydp, N_Vector x, void *user_data);
                int fdydx(realtype t, int it, realtype *dydx, N_Vector x, void *user_data);
                int fz(realtype t, int ie, int *nroots, realtype *z, N_Vector x, void *user_data);
                int fsz(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data);
                int fsz_tf(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data);
                int fdzdp(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data);
                int fdzdx(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data);
                int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data);
                int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data);
                int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data);
                int fdxdotdp(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data);
                int fdeltax(realtype t, int ie, realtype *deltax, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data);
                int fdeltasx(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data);
                int fdeltaxB(realtype t, int ie, realtype *deltaxB, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data);
                int fdeltaqB(realtype t, int ie, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data);
                int fsigma_y(realtype t, realtype *sigma_y, void *user_data);
                int fdsigma_ydp(realtype t, realtype *dsigma_ydp, void *user_data);
                int fsigma_z(realtype t, int ie, realtype *sigma_z, void *user_data);
                int fdsigma_zdp(realtype t, int ie, realtype *dsigma_zdp, void *user_data);
                int fJy(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data);
                int fJz(realtype t, int ie, realtype *Jz, realtype *z, N_Vector x, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data);
                int fdJydx(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data);
                int fdJydy(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data);
                int fdJydp(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data);
                int fdJzdx(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data);
                int fdJzdp(realtype t, int ie, realtype *dJzdp, realtype *z, N_Vector x, realtype *dzdp, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data);
                int fsJz(realtype t, int ie, realtype *sJz, realtype *s2Jz, realtype *dJzdz, realtype *dJzdp, realtype *sz, realtype *dzdp, realtype *mz, void *user_data, void *temp_data);
#endif /* _LW_cvodewrapfunctions */
