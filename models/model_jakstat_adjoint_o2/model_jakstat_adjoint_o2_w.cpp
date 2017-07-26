
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int w_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(udata->w,0,sizeof(realtype)*10);
  udata->w[0] = am_spline_pos(t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->w[1] = x_tmp[1]*x_tmp[1];
  udata->w[2] = 1.0/udata->k[0];
  udata->w[3] = 1.0/udata->k[1];
  udata->w[4] = x_tmp[3]*2.0;
  udata->w[5] = am_Dspline_pos(4,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->w[6] = am_Dspline_pos(6,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->w[7] = am_Dspline_pos(8,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->w[8] = am_Dspline_pos(10,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
  udata->w[9] = am_Dspline_pos(12,t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
return(status);

}


