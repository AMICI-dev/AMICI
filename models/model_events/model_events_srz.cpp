
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void srz_model_events(double *srz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip) {
switch (udata->plist[ip]) {
  case 0: {
    switch(ie) { 
        case 0: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 0)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 1)] = sx_tmp[0]-sx_tmp[2];

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
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 0)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 1)] = sx_tmp[0]-sx_tmp[2];

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
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 0)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 1)] = sx_tmp[0]-sx_tmp[2];

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
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 0)] = sx_tmp[1]-sx_tmp[2];

        } break;

        case 1: {
  rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + 1)] = sx_tmp[0]-sx_tmp[2];

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

}
}

