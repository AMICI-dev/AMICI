
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_dirac_w.h"

int dJydy_model_dirac(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydy,0,sizeof(realtype)*udata->ny*udata->nytrue*udata->nJ);
status = w_model_dirac(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydy[iy+(0)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
}
return(status);

}


