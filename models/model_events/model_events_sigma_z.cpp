
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void sigma_z_model_events(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
  tdata->sigmaz[0] = 1.0;
  tdata->sigmaz[1] = 1.0;
}

