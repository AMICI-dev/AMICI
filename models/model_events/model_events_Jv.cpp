
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void Jv_model_events(realtype t, N_Vector x, N_Vector xdot, N_Vector v, N_Vector Jv, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  Jv_tmp[0] = -tdata->h[3]*tdata->p[0]*v_tmp[0];
  Jv_tmp[1] = -tdata->p[2]*v_tmp[1]+tdata->p[1]*v_tmp[0]*exp(t*(-1.0/1.0E1));
  Jv_tmp[2] = -v_tmp[2];
}

