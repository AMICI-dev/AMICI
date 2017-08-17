
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include "model_dirac_w.h"

int deltasx_model_dirac(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
int ip;
memset(tdata->deltasx,0,sizeof(realtype)*2*udata->nplist);
status = w_model_dirac(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
              switch(ie) { 
              case 0: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 1: {
              switch(ie) { 
              case 0: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 2: {
              switch(ie) { 
              case 0: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 0: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  tdata->deltasx[ip*udata->nx + 0] = -udata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  tdata->deltasx[ip*udata->nx + 1] = -udata->stau[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

}
}
return(status);

}


