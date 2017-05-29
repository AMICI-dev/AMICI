
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int dsigma_ydp_model_dirac(realtype t, realtype *dsigma_ydp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int ip;
memset(dsigma_ydp,0,sizeof(realtype)*1*udata->nplist);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return(status);

}


