
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int sigma_y_model_jakstat_adjoint_o2(realtype t, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(sigma_y,0,sizeof(realtype)*54);
  sigma_y[0] = p[14];
  sigma_y[1] = p[15];
  sigma_y[2] = p[16];
  sigma_y[17] = 1.0;
  sigma_y[35] = 1.0;
  sigma_y[53] = 1.0;
return(status);

}


