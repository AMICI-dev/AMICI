
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_steadystate_w.h"

int sroot_model_steadystate(realtype t, int ie, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data, void *temp_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
TempData *tdata = (TempData*) temp_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
status = w_model_steadystate(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (udata->plist[ip]) {
}
}
return(status);

}


