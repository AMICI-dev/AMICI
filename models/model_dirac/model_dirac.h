#ifndef _am_model_dirac_h
#define _am_model_dirac_h

#include "model_dirac_J.h"
#include "model_dirac_JB.h"
#include "model_dirac_JBand.h"
#include "model_dirac_JBandB.h"
#include "model_dirac_JSparse.h"
#include "model_dirac_JSparseB.h"
#include "model_dirac_Jv.h"
#include "model_dirac_JvB.h"
#include "model_dirac_Jy.h"
#include "model_dirac_Jz.h"
#include "model_dirac_dJydp.h"
#include "model_dirac_dJydx.h"
#include "model_dirac_dJydy.h"
#include "model_dirac_dJzdp.h"
#include "model_dirac_dJzdx.h"
#include "model_dirac_deltaqB.h"
#include "model_dirac_deltasx.h"
#include "model_dirac_deltax.h"
#include "model_dirac_deltaxB.h"
#include "model_dirac_dsigma_ydp.h"
#include "model_dirac_dsigma_zdp.h"
#include "model_dirac_dwdp.h"
#include "model_dirac_dwdx.h"
#include "model_dirac_dxdotdp.h"
#include "model_dirac_dydp.h"
#include "model_dirac_dydx.h"
#include "model_dirac_dzdp.h"
#include "model_dirac_dzdx.h"
#include "model_dirac_qBdot.h"
#include "model_dirac_root.h"
#include "model_dirac_sJy.h"
#include "model_dirac_sJz.h"
#include "model_dirac_sigma_y.h"
#include "model_dirac_sigma_z.h"
#include "model_dirac_sroot.h"
#include "model_dirac_stau.h"
#include "model_dirac_sx0.h"
#include "model_dirac_sxdot.h"
#include "model_dirac_sz.h"
#include "model_dirac_sz_tf.h"
#include "model_dirac_w.h"
#include "model_dirac_x0.h"
#include "model_dirac_xBdot.h"
#include "model_dirac_xdot.h"
#include "model_dirac_y.h"
#include "model_dirac_z.h"

int J_model_dirac(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int JB_model_dirac(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int JBand_model_dirac(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int JBandB_model_dirac(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int JSparse_model_dirac(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int JSparseB_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int Jv_model_dirac(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);
int JvB_model_dirac(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB);
int Jy_model_dirac(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data);
int Jz_model_dirac(realtype t, int ie, realtype *Jz, realtype *z, N_Vector x, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data);
int dJydp_model_dirac(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data);
int dJydx_model_dirac(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data);
int dJydy_model_dirac(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data);
int dJzdp_model_dirac(realtype t, int ie, realtype *dJzdp, realtype *z, N_Vector x, realtype *dzdp, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data);
int dJzdx_model_dirac(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data);
int deltaqB_model_dirac(realtype t, int ie, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data);
int deltasx_model_dirac(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data);
int deltax_model_dirac(realtype t, int ie, realtype *deltax, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data);
int deltaxB_model_dirac(realtype t, int ie, realtype *deltaxB, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data);
int dsigma_ydp_model_dirac(realtype t, realtype *dsigma_ydp, void *user_data);
int dsigma_zdp_model_dirac(realtype t, int ie, realtype *dsigma_zdp, void *user_data);
int dwdp_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int dwdx_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int dxdotdp_model_dirac(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data);
int dydp_model_dirac(realtype t, int it, realtype *dydp, N_Vector x, void *user_data);
int dydx_model_dirac(realtype t, int it, realtype *dydx, N_Vector x, void *user_data);
int dzdp_model_dirac(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data);
int dzdx_model_dirac(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data);
int qBdot_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data);
int root_model_dirac(realtype t, N_Vector x, realtype *root, void *user_data);
int sJy_model_dirac(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y, realtype *sy, realtype *dydp, realtype *my, void *user_data);
int sJz_model_dirac(realtype t, int ie, realtype *sJz, realtype *s2Jz, realtype *dJzdz, realtype *dJzdp, realtype *sz, realtype *dzdp, realtype *mz, void *user_data, void *temp_data);
int sigma_y_model_dirac(realtype t, realtype *sigma_y, void *user_data);
int sigma_z_model_dirac(realtype t, int ie, realtype *sigma_z, void *user_data);
int sroot_model_dirac(realtype t, int ie, int *nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data);
int stau_model_dirac(realtype t, int ie, realtype *stau, N_Vector x, N_Vector *sx, void *user_data);
int sx0_model_dirac(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data);
int sxdot_model_dirac(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2);
int sz_model_dirac(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data);
int sz_tf_model_dirac(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data);
int w_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int x0_model_dirac(N_Vector x0, void *user_data);
int xBdot_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data);
int xdot_model_dirac(realtype t, N_Vector x, N_Vector xdot, void *user_data);
int y_model_dirac(realtype t, int it, realtype *y, N_Vector x, void *user_data);
int z_model_dirac(realtype t, int ie, int *nroots, realtype *z, N_Vector x, void *user_data);


#endif /* _am_model_dirac_h */
