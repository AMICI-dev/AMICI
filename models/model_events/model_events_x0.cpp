
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void x0_model_events(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0_tmp[0] = udata->k[0];
  x0_tmp[1] = udata->k[1];
  x0_tmp[2] = udata->k[2];
}

