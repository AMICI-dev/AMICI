
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_jakstat_adjoint_w.h"

int ddJydydy_model_jakstat_adjoint(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->ddJydydy,0,sizeof(realtype)*model->nytrue*model->ny*model->ny);
status = w_model_jakstat_adjoint(t,x,NULL,tdata);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->ddJydydy[iy+0*model->ny+0*model->ny*model->ny] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
}
if(!amiIsNaN(edata->my[1* udata->nt+it])){
    iy = 1;
  tdata->ddJydydy[iy+1*model->ny+1*model->ny*model->ny] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
}
if(!amiIsNaN(edata->my[2* udata->nt+it])){
    iy = 2;
  tdata->ddJydydy[iy+2*model->ny+2*model->ny*model->ny] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
}
return(status);

}


