
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include "model_events_w.h"

int stau_model_events(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
memset(tdata->stau,0,sizeof(realtype)*udata->nplist);
status = w_model_events(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

        } break;

        case 2: {

        } break;

        case 3: {
  tdata->stau[ip] = 1.0;

        } break;

    } 

  } break;

}
}
return(status);

}


