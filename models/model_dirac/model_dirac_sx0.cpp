
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_dirac_w.h"

using namespace amici;

int sx0_model_dirac(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
TempData *tdata = (TempData*) user_data;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *dx_tmp = nullptr;
if(dx)
    dx_tmp = N_VGetArrayPointer(dx);
realtype *sx0_tmp;
int ip;
realtype t = udata->tstart;
for(ip = 0; ip<udata->nplist; ip++) {
sx0_tmp = N_VGetArrayPointer(sx0[ip]);
memset(sx0_tmp,0,sizeof(realtype)*2);
switch (udata->plist[ip]) {
}
}
return(status);

}


