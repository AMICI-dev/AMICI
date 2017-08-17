
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include "model_dirac_w.h"

int stau_model_dirac(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
memset(udata->stau,0,sizeof(realtype)*udata->nplist);
status = w_model_dirac(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 1: {
    switch(ie) { 
        case 0: {
  udata->stau[ip] = 1.0;

        } break;

        case 1: {
  udata->stau[ip] = 1.0;

        } break;

    } 

  } break;

}
}
return(status);

}


