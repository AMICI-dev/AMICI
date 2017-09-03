
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_events_w.h"

int dJydp_model_events(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydp,0,sizeof(realtype)*udata->nytrue*udata->nplist*udata->nJ);
status = w_model_events(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0* udata->nt+it])){
    iy = 0;
  dJydp[iy+(0*4+3)*1] = dydp[3]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
}
return(status);

}


