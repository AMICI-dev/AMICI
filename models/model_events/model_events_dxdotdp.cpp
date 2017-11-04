
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dxdotdp_model_events(realtype t, N_Vector x,, void *user_data) {
switch (udata->plist[ip]) {
  case 0: {
  tdata->dxdotdp[0 + ip*model->nx] = -tdata->h[3]*x_tmp[0];

  } break;

  case 1: {
  tdata->dxdotdp[1 + ip*model->nx] = x_tmp[0]*exp(t*(-1.0/1.0E1));

  } break;

  case 2: {
  tdata->dxdotdp[1 + ip*model->nx] = -x_tmp[1];

  } break;

}
}

