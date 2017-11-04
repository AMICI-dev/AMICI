
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JSparseB_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
  JB->data[0] = tdata->h[3]*tdata->p[0];
  JB->data[1] = -tdata->p[1]*exp(t*(-1.0/1.0E1));
  JB->data[2] = tdata->p[2];
  JB->data[3] = 1.0;
}

