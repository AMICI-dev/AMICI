
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int dwdp_model_jakstat_adjoint(realtype t, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dwdp_tmp,0,sizeof(realtype)*5);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  dwdp_tmp[0] = am_Dspline_pos(4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[1] = am_Dspline_pos(6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[2] = am_Dspline_pos(8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[3] = am_Dspline_pos(10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  dwdp_tmp[4] = am_Dspline_pos(12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
return(status);

}


