
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int sx0_model_steadystate(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx0_tmp;
int ip;
realtype t = tstart;
for(ip = 0; ip<nplist; ip++) {
sx0_tmp = N_VGetArrayPointer(sx0[plist[ip]]);
memset(sx0_tmp,0,sizeof(realtype)*3);
switch (plist[ip]) {
}
}
return(status);

}


