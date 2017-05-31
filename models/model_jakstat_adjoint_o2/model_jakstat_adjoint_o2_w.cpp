
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int w_model_jakstat_adjoint_o2(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(w_tmp,0,sizeof(realtype)*10);
  w_tmp[0] = am_spline_pos(t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w_tmp[1] = x_tmp[1]*x_tmp[1];
  w_tmp[2] = 1.0/k[0];
  w_tmp[3] = 1.0/k[1];
  w_tmp[4] = x_tmp[3]*2.0;
  w_tmp[5] = am_Dspline_pos(4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w_tmp[6] = am_Dspline_pos(6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w_tmp[7] = am_Dspline_pos(8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w_tmp[8] = am_Dspline_pos(10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w_tmp[9] = am_Dspline_pos(12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
return(status);

}


