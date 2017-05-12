
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int deltaqB_model_jakstat_adjoint_o2(realtype t, int ie, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
realtype *xdot_old_tmp = N_VGetArrayPointer(xdot_old);
int ip;
memset(deltaqB,0,sizeof(realtype)*nplist*ng);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
}
}
return(status);

}


