
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int dydp_model_jakstat_adjoint(realtype t, int it, N_Vector x, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
  case 4: {
  tdata->dydp[ip*model->ny + 0] = -1.0/(tdata->p[4]*tdata->p[4])*tdata->p[13]*(x_tmp[1]+x_tmp[2]*2.0);
  tdata->dydp[ip*model->ny + 1] = -1.0/(tdata->p[4]*tdata->p[4])*tdata->p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0);

  } break;

  case 5: {
  tdata->dydp[ip*model->ny + 2] = am_Dspline_pos(4,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);

  } break;

  case 6: {
  tdata->dydp[ip*model->ny + 2] = am_Dspline_pos(6,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);

  } break;

  case 7: {
  tdata->dydp[ip*model->ny + 2] = am_Dspline_pos(8,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);

  } break;

  case 8: {
  tdata->dydp[ip*model->ny + 2] = am_Dspline_pos(10,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);

  } break;

  case 9: {
  tdata->dydp[ip*model->ny + 2] = am_Dspline_pos(12,t,5,0.0,tdata->p[5],5.0,tdata->p[6],1.0E1,tdata->p[7],2.0E1,tdata->p[8],6.0E1,tdata->p[9],0.0,0.0);

  } break;

  case 10: {
  tdata->dydp[ip*model->ny + 1] = 1.0;

  } break;

  case 11: {
  tdata->dydp[ip*model->ny + 0] = 1.0;

  } break;

  case 12: {
  tdata->dydp[ip*model->ny + 1] = (x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0)/tdata->p[4];

  } break;

  case 13: {
  tdata->dydp[ip*model->ny + 0] = (x_tmp[1]+x_tmp[2]*2.0)/tdata->p[4];

  } break;

}
}
return(status);

}


