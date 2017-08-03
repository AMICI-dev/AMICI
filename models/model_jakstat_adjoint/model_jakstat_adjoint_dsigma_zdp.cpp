
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_jakstat_adjoint_w.h"

int dsigma_zdp_model_jakstat_adjoint(realtype t, int ie, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
int ip;
memset(tdata->dsigmazdp,0,sizeof(realtype)*0*udata->nplist);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return(status);

}


