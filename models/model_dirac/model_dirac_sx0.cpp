
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int sx0_model_dirac(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx0_tmp;
int ip;
realtype t = udata->tstart;
for(ip = 0; ip<udata->nplist; ip++) {
sx0_tmp = N_VGetArrayPointer(sx0[ip]);
memset(sx0_tmp,0,sizeof(realtype)*2);
switch (udata->plist[ip]) {
}
}
return(status);

}


