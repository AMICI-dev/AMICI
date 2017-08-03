
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include "model_neuron_w.h"

int z_model_neuron(realtype t, int ie, N_Vector x, void *user_data, TempData *tdata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_neuron(t,x,NULL,user_data);
    switch(ie) { 
        case 0: {
  rdata->z[tdata->nroots[ie]+udata->nmaxevent*0] = t;

        } break;

    } 
return(status);

}


