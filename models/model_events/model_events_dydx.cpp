
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dydx_model_events(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  tdata->dydx[0+0*1] = tdata->p[3];
  tdata->dydx[0+1*1] = tdata->p[3];
  tdata->dydx[0+2*1] = tdata->p[3];
}

