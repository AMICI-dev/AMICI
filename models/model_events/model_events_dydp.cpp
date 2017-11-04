
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dydp_model_events(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) {
switch (udata->plist[ip]) {
  case 3: {
  tdata->dydp[ip*model->ny + 0] = x_tmp[0]+x_tmp[1]+x_tmp[2];

  } break;

}
}

