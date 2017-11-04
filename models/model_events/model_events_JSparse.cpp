
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JSparse_model_events(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  J->data[0] = -tdata->h[3]*tdata->p[0];
  J->data[1] = tdata->p[1]*exp(t*(-1.0/1.0E1));
  J->data[2] = -tdata->p[2];
  J->data[3] = -1.0;
}

