
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void dzdx_model_events(double *dzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  tdata->dzdx[0+1*2] = 1.0/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));
  tdata->dzdx[0+2*2] = -1.0/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));
  tdata->dzdx[1+0*2] = 1.0/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);
  tdata->dzdx[1+2*2] = -1.0/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);
}

