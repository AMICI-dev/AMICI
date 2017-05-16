
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int sigma_z_model_jakstat_adjoint_o2(realtype t, int ie, realtype *sigma_z, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
memset(sigma_z,0,sizeof(realtype)*0);
return(status);

}


