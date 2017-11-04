
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void sigma_y_model_events(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
  tdata->sigmay[0] = 1.0;
}

