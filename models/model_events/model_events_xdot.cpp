
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void xdot_model_events(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
  xdot_tmp[0] = -tdata->h[3]*tdata->p[0]*x_tmp[0];
  xdot_tmp[1] = -tdata->p[2]*x_tmp[1]+tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1));
  xdot_tmp[2] = tdata->h[2]-x_tmp[2];
}

