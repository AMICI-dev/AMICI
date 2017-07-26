
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int root_model_dirac(realtype t, N_Vector x, realtype *root, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_dirac(t,x,NULL,user_data);
  root[0] = -t+udata->p[1];
  root[1] = t-udata->p[1];
return(status);

}


