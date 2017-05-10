
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int y_model_jakstat_adjoint(realtype t, int it, realtype *y, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
  y[it+nt*(0)] = p[11]+(p[13]*(x_tmp[1]+x_tmp[2]*2.0))/p[4];
  y[it+nt*(1)] = p[10]+(p[12]*(x_tmp[0]+x_tmp[1]+x_tmp[2]*2.0))/p[4];
  y[it+nt*(2)] = am_spline_pos(t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
return(status);

}


