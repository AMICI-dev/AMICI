
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int x0_model_jakstat_adjoint_o2(N_Vector x0, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x0_tmp = N_VGetArrayPointer(x0);
memset(x0_tmp,0,sizeof(realtype)*162);
realtype t = udata->tstart;
  x0_tmp[0] = udata->p[4];
  x0_tmp[45] = 1.0;
return(status);

}


