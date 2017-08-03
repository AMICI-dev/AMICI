
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_neuron_o2_w.h"

int dJydsigma_model_neuron_o2(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydsigma,0,sizeof(realtype)*udata->nytrue*udata->ny*udata->nJ);
status = w_model_neuron_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydsigma[iy+(0+0*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*-1.0+1.0/tdata->sigmay[0];
  tdata->dJydsigma[iy+(1+0*5)*udata->nytrue] = rdata->y[it + udata->nt*1]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(2+0*5)*udata->nytrue] = rdata->y[it + udata->nt*2]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(3+0*5)*udata->nytrue] = rdata->y[it + udata->nt*3]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(4+0*5)*udata->nytrue] = rdata->y[it + udata->nt*4]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
}
return(status);

}


