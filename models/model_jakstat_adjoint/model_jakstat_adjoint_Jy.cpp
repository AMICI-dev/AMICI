
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int Jy_model_jakstat_adjoint(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  Jy[0] += amilog((sigma_y[0]*sigma_y[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[0]*sigma_y[0])*pow(my[it+nt*0]-y[it+nt*0],2.0)*5.0E-1;
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  Jy[0] += amilog((sigma_y[1]*sigma_y[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[1]*sigma_y[1])*pow(my[it+nt*1]-y[it+nt*1],2.0)*5.0E-1;
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  Jy[0] += amilog((sigma_y[2]*sigma_y[2])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[2]*sigma_y[2])*pow(my[it+nt*2]-y[it+nt*2],2.0)*5.0E-1;
}
return(status);

}


