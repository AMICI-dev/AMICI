
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int dzdp_model_dirac(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
int ip;
status = w_model_dirac(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
switch (plist[ip]) {
}
}
return(status);

}


