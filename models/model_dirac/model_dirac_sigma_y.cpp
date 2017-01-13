
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int sigma_y_model_dirac(realtype t, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(sigma_y,0,sizeof(realtype)*1);
  sigma_y[0] = 1.0;
return(status);

}


