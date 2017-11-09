
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include "model_neuron_o2_w.h"

using namespace amici;

void sz_model_neuron_o2(realtype t, int ie, N_Vector x, N_Vector *sx, amici::TempData *tdata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
w_model_neuron_o2(t,x,NULL,tdata);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 1)] = x_tmp[2]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))*-2.0-sx_tmp[2]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[2]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 2)] = -x_tmp[4]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[4]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[4]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 3)] = -x_tmp[6]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[6]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[6]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 4)] = -x_tmp[8]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[8]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[8]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 1)] = -x_tmp[4]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[2]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[2]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 2)] = x_tmp[4]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))*-2.0-sx_tmp[4]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[4]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 3)] = -x_tmp[6]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[4]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[6]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[6]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 4)] = -x_tmp[8]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[4]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[8]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[8]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 1)] = -x_tmp[6]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[2]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[2]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 2)] = -x_tmp[6]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[4]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[4]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[4]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 3)] = x_tmp[6]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))*-2.0-sx_tmp[6]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[6]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 4)] = -x_tmp[8]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[6]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[8]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[8]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 0)] = -sx_tmp[0]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 1)] = -x_tmp[8]*(x_tmp[3]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[2]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[2]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[2]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[2]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 2)] = -x_tmp[8]*(x_tmp[5]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[4]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[4]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[4]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[4]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 3)] = -x_tmp[8]*(x_tmp[7]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-x_tmp[6]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))-sx_tmp[6]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[6]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[6]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);
  rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*model->nz + 4)] = x_tmp[8]*(x_tmp[9]*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)-x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0)*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0))*-2.0-sx_tmp[8]/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)-(sx_tmp[0]*((x_tmp[8]*(x_tmp[0]*(2.0/2.5E1)+5.0))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2)+tdata->p[0]*x_tmp[8]*(x_tmp[1]-tdata->p[1]*x_tmp[0])*1.0/pow(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2,2.0)))/(udata->k[1]+x_tmp[0]*5.0-x_tmp[1]+(x_tmp[0]*x_tmp[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

}
}
return;

}


