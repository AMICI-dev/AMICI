
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_w.h"

int Jv_model_dirac(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xdot_tmp = N_VGetArrayPointer(xdot);
realtype *v_tmp = N_VGetArrayPointer(v);
realtype *Jv_tmp = N_VGetArrayPointer(Jv);
memset(Jv_tmp,0,sizeof(realtype)*2);
status = w_model_dirac(t,x,NULL,user_data);
  Jv_tmp[0] = -udata->p[0]*v_tmp[0];
  Jv_tmp[1] = udata->p[2]*v_tmp[0]-udata->p[3]*v_tmp[1];
return(status);

}


