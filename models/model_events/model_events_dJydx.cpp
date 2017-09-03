
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_events_w.h"

int dJydx_model_events(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_events(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0* udata->nt+it])){
    iy = 0;
  dJydx[it+(0+0*3)*udata->nt] += dydx[0]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydx[it+(1+0*3)*udata->nt] += dydx[1]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydx[it+(2+0*3)*udata->nt] += dydx[2]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
}
return(status);

}


