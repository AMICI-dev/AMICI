
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int dJydy_model_dirac(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydy,0,sizeof(realtype)*nytrue*nytrue*ng);
status = w_model_dirac(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  dJydy[0] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
}
return(status);

}


