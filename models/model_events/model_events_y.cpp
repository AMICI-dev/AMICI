
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void y_model_events(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  rdata->y[it + udata->nt*0] = tdata->p[3]*(x_tmp[0]+x_tmp[1]+x_tmp[2]);
}

