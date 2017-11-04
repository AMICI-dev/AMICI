
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void root_model_events(realtype t, N_Vector x, realtype *root, void *user_data) {
  root[0] = x_tmp[1]-x_tmp[2];
  root[1] = x_tmp[0]-x_tmp[2];
  root[2] = t-4.0;
  root[3] = t-tdata->p[3];
}

