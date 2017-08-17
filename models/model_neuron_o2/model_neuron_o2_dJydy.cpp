
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_neuron_o2_w.h"

int dJydy_model_neuron_o2(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydy,0,sizeof(realtype)*udata->ny*udata->nytrue*udata->nJ);
status = w_model_neuron_o2(t,x,NULL,tdata);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydy[iy+(0+0*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(1+0*5)*udata->nytrue] = rdata->y[it + udata->nt*1]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(1+1*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(2+0*5)*udata->nytrue] = rdata->y[it + udata->nt*2]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(2+2*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(3+0*5)*udata->nytrue] = rdata->y[it + udata->nt*3]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(3+3*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(4+0*5)*udata->nytrue] = rdata->y[it + udata->nt*4]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(4+4*5)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
}
return(status);

}


