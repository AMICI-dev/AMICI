
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_jakstat_adjoint_w.h"

int y_model_jakstat_adjoint(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  rdata->y[it + udata->nt*0] = udata->p[11]+(udata->p[13]*(x_tmp[1]+x_tmp[2]*2.0))/udata->p[4];
  rdata->y[it + udata->nt*1] = udata->p[10]+(udata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/udata->p[4];
  rdata->y[it + udata->nt*2] = am_spline_pos(t,5,0.0,udata->p[5],5.0,udata->p[6],1.0E1,udata->p[7],2.0E1,udata->p[8],6.0E1,udata->p[9],0.0,0.0);
return(status);

}


