
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include <include/tdata.h>
#include <include/tdata_accessors.h>
#undef t
#undef x
#undef x_tmp
#undef dzdp
#undef dzdx
#undef dx
#undef sigma_y
#undef sigma_z
#undef dsigma_ydp
#undef dsigma_zdp
#include "model_jakstat_adjoint_o2_w.h"

int dJzdx_model_jakstat_adjoint_o2(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
return(status);

}


