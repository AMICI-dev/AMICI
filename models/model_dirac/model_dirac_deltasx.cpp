
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int deltasx_model_dirac(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
int ip;
memset(deltasx,0,sizeof(realtype)*2*nplist);
status = w_model_dirac(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (plist[ip]) {
  case 0: {
              switch(ie) { 
              case 0: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 1: {
              switch(ie) { 
              case 0: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 2: {
              switch(ie) { 
              case 0: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 0: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              case 1: {
  deltasx[ip*2 + 0] = -stau_tmp[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);
  deltasx[ip*2 + 1] = -stau_tmp[ip]*(xdot_tmp[1]-xdot_old_tmp[1]);

              } break;

              } 

  } break;

}
}
return(status);

}


