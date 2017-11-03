
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_nested_events_w.h"

using namespace amici;

void stau_model_nested_events(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
memset(tdata->stau,0,sizeof(realtype)*udata->nplist);
w_model_nested_events(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {

        } break;

        case 1: {
  tdata->stau[ip] = sx_tmp[0]/(tdata->p[4]*x_tmp[0]-tdata->h[1]*tdata->p[3]*x_tmp[0]);

        } break;

        case 2: {

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {

        } break;

        case 1: {
  tdata->stau[ip] = sx_tmp[0]/(tdata->p[4]*x_tmp[0]-tdata->h[1]*tdata->p[3]*x_tmp[0]);

        } break;

        case 2: {

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = 1.0;

        } break;

        case 1: {
  tdata->stau[ip] = sx_tmp[0]/(tdata->p[4]*x_tmp[0]-tdata->h[1]*tdata->p[3]*x_tmp[0]);

        } break;

        case 2: {
  tdata->stau[ip] = 1.0;

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {

        } break;

        case 1: {
  tdata->stau[ip] = sx_tmp[0]/(tdata->p[4]*x_tmp[0]-tdata->h[1]*tdata->p[3]*x_tmp[0]);

        } break;

        case 2: {

        } break;

    } 

  } break;

  case 4: {
    switch(ie) { 
        case 0: {

        } break;

        case 1: {
  tdata->stau[ip] = sx_tmp[0]/(tdata->p[4]*x_tmp[0]-tdata->h[1]*tdata->p[3]*x_tmp[0]);

        } break;

        case 2: {

        } break;

    } 

  } break;

}
}
return;

}


