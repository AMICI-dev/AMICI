#ifndef _am_model_dirac_h
#define _am_model_dirac_h

#include "model_dirac_J.h"
#include "model_dirac_JB.h"
#include "model_dirac_JBand.h"
#include "model_dirac_JBandB.h"
#include "model_dirac_JDiag.h"
#include "model_dirac_JSparse.h"
#include "model_dirac_JSparseB.h"
#include "model_dirac_Jrz.h"
#include "model_dirac_Jv.h"
#include "model_dirac_JvB.h"
#include "model_dirac_Jy.h"
#include "model_dirac_Jz.h"
#include "model_dirac_dJrzdsigma.h"
#include "model_dirac_dJrzdz.h"
#include "model_dirac_dJydsigma.h"
#include "model_dirac_dJydy.h"
#include "model_dirac_dJzdsigma.h"
#include "model_dirac_dJzdz.h"
#include "model_dirac_deltaqB.h"
#include "model_dirac_deltasx.h"
#include "model_dirac_deltax.h"
#include "model_dirac_deltaxB.h"
#include "model_dirac_drzdp.h"
#include "model_dirac_drzdx.h"
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
#include "model_dirac_rz.h"
#include "model_dirac_sigma_y.h"
#include "model_dirac_sigma_z.h"
#include "model_dirac_srz.h"
#include "model_dirac_stau.h"
#include "model_dirac_sx0.h"
#include "model_dirac_sxdot.h"
#include "model_dirac_sz.h"
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
int JDiag_model_dirac(realtype t, N_Vector JDiag, N_Vector x, void *user_data);
int JSparse_model_dirac(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int JSparseB_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
int Jrz_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int Jv_model_dirac(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);
int JvB_model_dirac(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB);
int Jy_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int Jz_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJrzdsigma_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJrzdz_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJydsigma_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJydy_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJzdsigma_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int dJzdz_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata);
int deltaqB_model_dirac(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata);
int deltasx_model_dirac(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data, TempData *tdata);
int deltax_model_dirac(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata);
int deltaxB_model_dirac(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data, TempData *tdata);
int drzdp_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata);
int drzdx_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata);
int dsigma_ydp_model_dirac(realtype t, void *user_data, TempData *tdata);
int dsigma_zdp_model_dirac(realtype t, int ie, void *user_data, TempData *tdata);
int dwdp_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int dwdx_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int dxdotdp_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int dydp_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata);
int dydx_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata);
int dzdp_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata);
int dzdx_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata);
int qBdot_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data);
int root_model_dirac(realtype t, N_Vector x, realtype *root, void *user_data);
int rz_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata);
int sigma_y_model_dirac(realtype t, void *user_data, TempData *tdata);
int sigma_z_model_dirac(realtype t, int ie, void *user_data, TempData *tdata);
int srz_model_dirac(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata);
int stau_model_dirac(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata);
int sx0_model_dirac(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data);
int sxdot_model_dirac(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2);
int sz_model_dirac(realtype t, int ie, N_Vector x, N_Vector *sx, void *user_data, TempData *tdata, ReturnData *rdata);
int w_model_dirac(realtype t, N_Vector x, N_Vector dx, void *user_data);
int x0_model_dirac(N_Vector x0, void *user_data);
int xBdot_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data);
int xdot_model_dirac(realtype t, N_Vector x, N_Vector xdot, void *user_data);
int y_model_dirac(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata);
int z_model_dirac(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata);


#endif /* _am_model_dirac_h */
