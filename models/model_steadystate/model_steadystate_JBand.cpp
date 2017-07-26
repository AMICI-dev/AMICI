
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_steadystate_J.h"
#include "model_steadystate_w.h"

int JBand_model_steadystate(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
int status = 0;
return(J_model_steadystate(N, t, x, xdot, J, user_data, tmp1, tmp2, tmp3));}


