
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void xBdot_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) {
  xBdot_tmp[0] = -tdata->p[1]*xB_tmp[1]*exp(t*(-1.0/1.0E1))+tdata->h[3]*tdata->p[0]*xB_tmp[0];
  xBdot_tmp[1] = tdata->p[2]*xB_tmp[1];
  xBdot_tmp[2] = xB_tmp[2];
}

