
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void z_model_events(double *z, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
    switch(ie) { 
        case 0: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*0] = t;

        } break;

        case 1: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*1] = t;

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 
}

