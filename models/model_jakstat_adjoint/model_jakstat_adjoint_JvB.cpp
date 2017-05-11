
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int JvB_model_jakstat_adjoint(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*9);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  JvB_tmp[0] = p[0]*vB_tmp[0]*w_tmp[0]-p[0]*vB_tmp[1]*w_tmp[0];
  JvB_tmp[1] = dwdx_tmp[0]*p[1]*vB_tmp[1]*2.0-dwdx_tmp[0]*p[1]*vB_tmp[2];
  JvB_tmp[2] = p[2]*vB_tmp[2]-(k[0]*p[2]*vB_tmp[3])/k[1];
  JvB_tmp[3] = p[3]*vB_tmp[3]-p[3]*vB_tmp[4]*2.0;
  JvB_tmp[4] = p[3]*vB_tmp[4]-p[3]*vB_tmp[5];
  JvB_tmp[5] = p[3]*vB_tmp[5]-p[3]*vB_tmp[6];
  JvB_tmp[6] = p[3]*vB_tmp[6]-p[3]*vB_tmp[7];
  JvB_tmp[7] = p[3]*vB_tmp[7]-p[3]*vB_tmp[8];
  JvB_tmp[8] = p[3]*vB_tmp[8]-(k[1]*p[3]*vB_tmp[0])/k[0];
return(status);

}


