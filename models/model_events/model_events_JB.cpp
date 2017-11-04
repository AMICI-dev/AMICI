
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JB_model_events(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
  JB->data[0+0*3] = tdata->h[3]*tdata->p[0];
  JB->data[0+1*3] = -tdata->p[1]*exp(t*(-1.0/1.0E1));
  JB->data[1+1*3] = tdata->p[2];
  JB->data[2+2*3] = 1.0;
}

