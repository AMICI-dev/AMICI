
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

using namespace amici;

void dzdp_model_events(realtype t, int ie, N_Vector x, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
int ip;
w_model_events(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return;

}


