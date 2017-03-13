
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int deltax_model_dirac(realtype t, int ie, realtype *deltax, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
memset(deltax,0,sizeof(realtype)*2);
status = w_model_dirac(t,x,NULL,user_data);
              switch(ie) { 
              case 0: {
  deltax[0] = 1.0;

              } break;

              case 1: {
  deltax[0] = 1.0;
  deltax[1] = 1.0;

              } break;

              } 
return(status);

}


