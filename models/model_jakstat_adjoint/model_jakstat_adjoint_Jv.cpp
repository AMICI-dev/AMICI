
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int Jv_model_jakstat_adjoint(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*9);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  Jv_tmp[0] = -p[0]*v_tmp[0]*w_tmp[0]+(k[1]*p[3]*v_tmp[8])/k[0];
  Jv_tmp[1] = dwdx_tmp[0]*p[1]*v_tmp[1]*-2.0+p[0]*v_tmp[0]*w_tmp[0];
  Jv_tmp[2] = -p[2]*v_tmp[2]+dwdx_tmp[0]*p[1]*v_tmp[1];
  Jv_tmp[3] = -p[3]*v_tmp[3]+(k[0]*p[2]*v_tmp[2])/k[1];
  Jv_tmp[4] = p[3]*v_tmp[3]*2.0-p[3]*v_tmp[4];
  Jv_tmp[5] = p[3]*v_tmp[4]-p[3]*v_tmp[5];
  Jv_tmp[6] = p[3]*v_tmp[5]-p[3]*v_tmp[6];
  Jv_tmp[7] = p[3]*v_tmp[6]-p[3]*v_tmp[7];
  Jv_tmp[8] = p[3]*v_tmp[7]-p[3]*v_tmp[8];
return(status);

}


