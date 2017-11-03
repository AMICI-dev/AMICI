
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_dirac_w.h"

using namespace amici;

void dJydsigma_model_dirac(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydsigma,0,sizeof(realtype)*model->nytrue*model->ny*model->nJ);
w_model_dirac(t,x,NULL,tdata);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydsigma[iy+(0)*model->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0]*1.0,2.0)*-1.0+1.0/tdata->sigmay[0];
}
return;

}


