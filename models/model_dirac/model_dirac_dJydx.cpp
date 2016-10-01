
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int dJydx_model_dirac(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_dirac(t,x,NULL,user_data);
if(!amiIsNaN(my[0*nt+it])){
  dJydx[it+(1+0*2)*nt] += dydx[1]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
}
return(status);

}


