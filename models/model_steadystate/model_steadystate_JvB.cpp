
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_w.h"

int JvB_model_steadystate(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
realtype *vB_tmp = N_VGetArrayPointer(vB);
realtype *JvB_tmp = N_VGetArrayPointer(JvB);
memset(JvB_tmp,0,sizeof(realtype)*3);
status = w_model_steadystate(t,x,NULL,user_data);
  JvB_tmp[0] = vB_tmp[0]*(udata->p[1]*x_tmp[1]+udata->p[0]*udata->dwdx[0]*2.0)+vB_tmp[1]*(udata->p[1]*x_tmp[1]-udata->p[0]*udata->dwdx[0])-udata->p[1]*x_tmp[1]*vB_tmp[2];
  JvB_tmp[1] = -vB_tmp[0]*(udata->p[2]*2.0-udata->p[1]*x_tmp[0])+vB_tmp[1]*(udata->p[2]+udata->p[1]*x_tmp[0])-udata->p[1]*x_tmp[0]*vB_tmp[2];
  JvB_tmp[2] = vB_tmp[2]*(udata->k[3]+udata->dwdx[1])-vB_tmp[0]*udata->dwdx[1]-vB_tmp[1]*udata->dwdx[1];
return(status);

}


