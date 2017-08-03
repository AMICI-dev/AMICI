
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_dwdx.h"
#include "model_neuron_w.h"

int JSparseB_model_neuron(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
realtype *xB_tmp = N_VGetArrayPointer(xB);
realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);
  SparseSetMatToZero(JB);
  JB->indexvals[0] = 0;
  JB->indexvals[1] = 1;
  JB->indexvals[2] = 0;
  JB->indexvals[3] = 1;
  JB->indexptrs[0] = 0;
  JB->indexptrs[1] = 2;
  JB->indexptrs[2] = 4;
status = w_model_neuron(t,x,NULL,user_data);
status = dwdx_model_neuron(t,x,NULL,user_data);
  JB->data[0] = x_tmp[0]*(-2.0/2.5E1)-5.0;
  JB->data[1] = 1.0;
  JB->data[2] = -udata->p[0]*udata->p[1];
  JB->data[3] = udata->p[0];
return(status);

}


