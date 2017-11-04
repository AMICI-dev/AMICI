
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JvB_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, N_Vector vB, N_Vector JvB, void *user_data, N_Vector tmpB1, N_Vector tmpB2) {
  JvB_tmp[0] = -tdata->p[1]*vB_tmp[1]*exp(t*(-1.0/1.0E1))+tdata->h[3]*tdata->p[0]*vB_tmp[0];
  JvB_tmp[1] = tdata->p[2]*vB_tmp[1];
  JvB_tmp[2] = vB_tmp[2];
}

