
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int dJydx_model_steadystate(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_steadystate(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  dJydx[it+(0+0*3)*nt] += dydx[0]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  dJydx[it+(1+0*3)*nt] += dydx[4]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  dJydx[it+(2+0*3)*nt] += dydx[8]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
}
return(status);

}


