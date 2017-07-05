
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_steadystate_w.h"

int dJzdp_model_steadystate(realtype t, int ie, N_Vector x, realtype *z, realtype *mz, realtype *dzdp, void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_steadystate(t,x,NULL,user_data);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return(status);

}


