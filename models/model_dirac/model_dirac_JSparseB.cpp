
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_dirac_dwdx.h"
#include "model_dirac_w.h"

int JSparseB_model_dirac(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  SparseSetMatToZero(JB);
  JB->indexvals[0] = 0;
  JB->indexvals[1] = 0;
  JB->indexvals[2] = 1;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 1;
  JB->indexptrs[2] = 3;
status = w_model_dirac(t,x,NULL,user_data);
status = dwdx_model_dirac(t,x,NULL,user_data);
  JB->data[0] = udata->p[0];
  JB->data[1] = -udata->p[2];
  JB->data[2] = udata->p[3];
return(status);

}


