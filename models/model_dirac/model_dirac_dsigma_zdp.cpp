
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int dsigma_zdp_model_dirac(realtype t, int ie, realtype *dsigma_zdp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int ip;
memset(dsigma_zdp,0,sizeof(realtype)*0*udata->nplist);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return(status);

}


