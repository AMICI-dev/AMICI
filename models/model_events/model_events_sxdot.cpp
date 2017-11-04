
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void sxdot_model_events(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  sxdot_tmp[0] = tdata->dxdotdp[0 + ip*model->nx]+tdata->J->data[0]*sx_tmp[0];
  sxdot_tmp[1] = tdata->dxdotdp[1 + ip*model->nx]+tdata->J->data[1]*sx_tmp[0]+tdata->J->data[2]*sx_tmp[1];
  sxdot_tmp[2] = tdata->dxdotdp[2 + ip*model->nx]+tdata->J->data[3]*sx_tmp[2];
}

