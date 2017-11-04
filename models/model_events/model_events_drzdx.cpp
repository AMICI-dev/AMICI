
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void drzdx_model_events(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  tdata->drzdx[0+0*1] = 1.0;
  tdata->drzdx[0+2*1] = -1.0;
}

