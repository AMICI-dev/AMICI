
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

int dJzdsigma_model_jakstat_adjoint_o2(realtype t, int ie, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJzdsigma,0,sizeof(realtype)*model->nztrue*model->nz*model->nJ);
status = w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
return(status);

}


