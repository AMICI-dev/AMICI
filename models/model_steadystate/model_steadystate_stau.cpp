
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int stau_model_steadystate(realtype t, int ie, realtype *stau, N_Vector x, N_Vector *sx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *sx_tmp;
int ip;
memset(stau,0,sizeof(realtype)*nplist);
status = w_model_steadystate(t,x,NULL,user_data);
for(ip = 0; ip<nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (plist[ip]) {
}
}
return(status);

}


