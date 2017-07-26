
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_neuron_o2_JB.h"
#include "model_neuron_o2_w.h"

int JBandB_model_neuron_o2(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {
int status = 0;
return(JB_model_neuron_o2(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B));}


