
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

using namespace amici;

int deltasx_model_events(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, amici::TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
realtype *xdot_tmp = nullptr;
if(xdot)
    xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = nullptr;
if(xdot_old)
    xdot_old_tmp = N_VGetArrayPointer(xdot_old);
int ip;
memset(tdata->deltasx,0,sizeof(realtype)*3*udata->nplist);
status = w_model_events(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 1: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 2: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

}
}
return(status);

}


