
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int sx0_model_dirac(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx0_tmp;
int ip;
for(ip = 0; ip<np; ip++) {
sx0_tmp = N_VGetArrayPointer(sx0[plist[ip]]);
memset(sx0_tmp,0,sizeof(realtype)*2);
switch (plist[ip]) {
}
}
return(status);

}


