
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include "model_neuron_w.h"

int dsigma_ydp_model_neuron(realtype t, void *user_data, TempData *tdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
int ip;
memset(tdata->dsigmaydp,0,sizeof(realtype)*1*udata->nplist);
for(ip = 0; ip<udata->nplist; ip++) {
switch (udata->plist[ip]) {
}
}
return(status);

}


