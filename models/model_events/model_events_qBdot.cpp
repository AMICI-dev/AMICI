
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void qBdot_model_events(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) {
switch (udata->plist[ip]) {
  case 0: {
  qBdot_tmp[ip + udata->nplist*0] = tdata->h[3]*x_tmp[0]*xB_tmp[0];

  } break;

  case 1: {
  qBdot_tmp[ip + udata->nplist*0] = -x_tmp[0]*xB_tmp[1]*exp(t*(-1.0/1.0E1));

  } break;

  case 2: {
  qBdot_tmp[ip + udata->nplist*0] = x_tmp[1]*xB_tmp[1];

  } break;

}
}

