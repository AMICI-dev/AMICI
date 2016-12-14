
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int Jy_model_dirac(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_dirac(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  Jy[0] += amilog((sigma_y[0]*sigma_y[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[0]*sigma_y[0])*pow(my[it+nt*0]-y[it+nt*0],2.0)*5.0E-1;
}
return(status);

}


