
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int sigma_y_model_jakstat_adjoint(realtype t, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(sigma_y,0,sizeof(realtype)*3);
  sigma_y[0] = p[14];
  sigma_y[1] = p[15];
  sigma_y[2] = p[16];
return(status);

}


