
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void JDiag_model_events(realtype t, N_Vector JDiag, N_Vector x, void *user_data) {
  JDiag_tmp[0+0*3] = -tdata->h[3]*tdata->p[0];
  JDiag_tmp[1+0*3] = -tdata->p[2];
  JDiag_tmp[2+0*3] = -1.0;
}

