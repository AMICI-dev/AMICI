
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
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
  JvB_tmp[0] = udata->p[0]*udata->w[0]*vB_tmp[0]-udata->p[0]*udata->w[0]*vB_tmp[1];
  JvB_tmp[1] = udata->p[1]*vB_tmp[1]*udata->dwdx[0]*2.0-udata->p[1]*vB_tmp[2]*udata->dwdx[0];
  JvB_tmp[2] = udata->p[2]*vB_tmp[2]-(udata->k[0]*udata->p[2]*vB_tmp[3])/udata->k[1];
  JvB_tmp[3] = udata->p[3]*vB_tmp[3]-udata->p[3]*vB_tmp[4]*2.0;
  JvB_tmp[4] = udata->p[3]*vB_tmp[4]-udata->p[3]*vB_tmp[5];
  JvB_tmp[5] = udata->p[3]*vB_tmp[5]-udata->p[3]*vB_tmp[6];
  JvB_tmp[6] = udata->p[3]*vB_tmp[6]-udata->p[3]*vB_tmp[7];
  JvB_tmp[7] = udata->p[3]*vB_tmp[7]-udata->p[3]*vB_tmp[8];
  JvB_tmp[8] = udata->p[3]*vB_tmp[8]-(udata->k[1]*udata->p[3]*vB_tmp[0])/udata->k[0];
return(status);

}


