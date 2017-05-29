
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int deltaxB_model_jakstat_adjoint_o2(realtype t, int ie, realtype *deltaxB, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
memset(deltaxB,0,sizeof(realtype)*162);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
return(status);

}


