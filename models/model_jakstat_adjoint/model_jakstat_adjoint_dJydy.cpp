
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int dJydy_model_jakstat_adjoint(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydy,0,sizeof(realtype)*nytrue*nytrue*ng);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  dJydy[iy+(0*1+0)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  dJydy[iy+(1*1+0)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  dJydy[iy+(2*1+0)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
}
return(status);

}


