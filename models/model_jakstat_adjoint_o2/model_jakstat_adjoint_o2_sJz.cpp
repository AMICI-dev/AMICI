
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

int sJz_model_jakstat_adjoint_o2(realtype t, int ie, realtype *sJz, realtype *s2Jz, realtype *dJzdz, realtype *dJzdp, realtype *sz, realtype *dzdp, realtype *mz, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
int ip;
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
}
}
return(status);

}


